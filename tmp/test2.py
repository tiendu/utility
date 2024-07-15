import argparse
import os
import gzip
import logging
import hashlib
import random
from collections import Counter
from functools import partial
from itertools import groupby
from dataclasses import dataclass
from typing import List
import concurrent.futures

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna']
LINE_LENGTH = 80
CHUNK_SIZE = 100_000

# Setup logging to display messages with INFO level and above.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Data class to represent a sequence.
@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __hash__(self):
        '''Define a custom hash function for the dataclass.'''
        return hash((self.id, self.sequence, self.quality))

def read_sequences_from_file(file_path: str, file_type: str) -> List[Seq]:
    '''Read sequences from a file and return them as a list of Seq objects.'''
    sequences = []

    # Determine the appropriate file opening mode based on the file extension.
    opener = gzip.open if file_path.endswith('.gz') else open

    if file_type == 'FASTQ':
        with opener(file_path, 'rt') as fin:
            # Group lines into sets of 4 (header, sequence, '+', quality).
            groups = groupby(enumerate(fin), key=lambda x: x[0] // 4)
            for _, group in groups:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                name = header_line[1:]  # remove '@' character
                seq = sequence_line.upper()
                qual = quality_line.upper()
                sequences.append(Seq(name, seq, qual))
    elif file_type == 'FASTA':
        with opener(file_path, 'rt') as fin:
            # Group lines into sequences based on '>' character indicating header lines.
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == '>'))
            for header in faiter:
                headerstr = next(header).strip()
                name = headerstr[1:]  # remove '>' character
                seq = ''.join(s.strip().upper() for s in next(faiter))
                sequences.append(Seq(name, seq))
    return sequences

def write_sequences_to_file(sequences: List[Seq], file_path: str) -> None:
    '''Write sequences to a file.'''
    opener = gzip.open if file_path.endswith('.gz') else open
    
    with opener(file_path, 'wt') as f:
        for seq in sequences:
            if seq.quality == '':
                f.write(f'>{seq.id}\n')
                # Write sequence with a defined line length
                for i in range(0, len(seq.sequence), LINE_LENGTH):
                    f.write(seq.sequence[i:i+LINE_LENGTH] + '\n')
            else:
                f.write(f'@{seq.id}\n{seq.sequence}\n+\n{seq.quality}\n')

def hash_sequence(sequence: str, hash_function=hashlib.sha3_256) -> str:
    '''Return the hash of a sequence using the specified hash function.'''
    return hash_function(sequence.encode()).hexdigest()

def generate_kmers(string: str, k: int) -> List[str]:
    '''Split a string into k-mers of length k.'''
    return [string[i:i+k] for i in range(len(string) - k + 1)]

def deduplicate_chunk(sequences: List[Seq], uniq_seqs: dict) -> List[Seq]:
    '''Deduplicate sequences within a chunk.'''
    logging.info(f'Deduplicating a chunk with {len(sequences)} sequences...')
    sequences.sort(key=lambda s: len(s.sequence), reverse=True)
    uniq_kmer_hashes = set()

    if sequences:
        min_length = len(sequences[-1].sequence)
    logging.info(f'Using k-mers of minimum length: {min_length}')

    for current_seq in sequences:
        kmers = generate_kmers(current_seq.sequence, min_length)
        kmer_hashes = {hash_sequence(kmer) for kmer in kmers}
        if all(hash_val in uniq_kmer_hashes for hash_val in kmer_hashes):
            continue

        uniq_seqs[hash_sequence(current_seq.sequence)] = current_seq
        uniq_kmer_hashes.update(kmer_hashes)

    logging.info(f'Chunk deduplication complete. Unique sequences: {len(uniq_seqs)}')

    return list(uniq_seqs.values())
        
def deduplicate_concurrently(sequences: List[Seq], num_threads: int) -> List[Seq]:
    '''Perform recursive deduplication of sequences using multiple threads.'''
    while True:
        if not sequences:
            break

        logging.info(f'Initial number of sequences: {len(sequences)}')

        total_sequences = len(sequences)
        chunk_size = max(1, min(CHUNK_SIZE, total_sequences // num_threads))
        shared_sequences = dict()
        
        # Process each chunk concurrently using multiple threads.
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            func = partial(deduplicate_chunk, uniq_seqs=shared_sequences)
            futures = [executor.submit(func, sequences[i:i + chunk_size]) for i in range(0, total_sequences, chunk_size)]
            concurrent.futures.wait(futures)
            for future in concurrent.futures.as_completed(futures):
                for sequence in future.result():
                    shared_sequences[hash_sequence(sequence.sequence)] = sequence

        # If no sequences were deduplicated in this iteration, we're done.
        if len(shared_sequences) == len(sequences):
            break

        # Otherwise, shuffle the partially deduplicated sequences and repeat the process.
        sequences = list(shared_sequences.values())
        random.shuffle(sequences)
    
    return list(shared_sequences.values())

def main():
    '''Main function that handles command-line arguments and invokes the deduplication.'''
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')

    args = parser.parse_args()

    # Validate input and output file paths.
    if not os.path.exists(args.input_file):
        raise ValueError(f'The input file "{args.input_file}" does not exist.')
    try:
        with open(args.output_file, 'w'):
            pass
    except Exception as e:
        raise ValueError(f'Error creating output file "{args.output_file}": {e}')

    # Check the file type.
    file_type = ''
    if any(ext in args.input_file for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(ext in args.input_file for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).')

    # Validate the number of threads.
    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f'Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}')
        args.num_threads = min(max(args.num_threads, 1), max_threads)

    # Read sequences from the input file.
    logging.info(f'Reading sequences from file: {args.input_file}...')
    sequences = read_sequences_from_file(args.input_file, file_type)

    # Re-adjust the number of threads according to the number of sequences.
    if args.num_threads >= len(sequences) * 0.1:
        logging.warning(f'Number of sequences too low, thread count adjusted to 1')
        args.num_threads = 1

    # Perform deduplication.
    deduplicated_sequences = deduplicate_concurrently(sequences, args.num_threads)

    # Write deduplicated sequences to the output file.
    logging.info(f'Finished. Wrote {len(deduplicated_sequences)} final deduplicated sequences to: {args.output_file}')
    write_sequences_to_file(deduplicated_sequences, args.output_file)

if __name__ == '__main__':
    main()
