import argparse
import os
import re
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

# Setup logging to display messages with INFO level and above.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Data class to represent a sequence.
@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __hash__(self):
        return hash((self.id, self.sequence, self.quality))

def read_sequences_from_file(file_path: str, file_type: str) -> List[Seq]:
    sequences = []

    # Determine the appropriate file opening mode based on the file extension.
    opener = gzip.open if file_path.endswith('.gz') else open

    if file_type == 'FASTQ':
        with opener(file_path, 'rt') as fin:
            groups = groupby(enumerate(fin), key=lambda x: x[0] // 4)
            for _, group in groups:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                name = header_line[1:]  # remove '@' character
                seq = sequence_line.upper()
                qual = quality_line.upper()
                sequences.append(Seq(name, seq, qual))
    elif file_type == 'FASTA':
        with opener(file_path, 'rt') as fin:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == '>'))
            for header in faiter:
                headerstr = next(header).strip()
                name = headerstr[1:]  # remove '>' character
                seq = ''.join(s.strip().upper() for s in next(faiter))
                sequences.append(Seq(name, seq))
    return sequences

def write_sequences_to_file(sequences: List[Seq], file_path: str) -> None:
    with open(file_path, 'w') as f:
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

def deduplicate_chunk(sequences: List[Seq], unique_seqs: dict) -> List[Seq]:
    logging.info(f'Processing a chunk with {len(sequences)} sequences')
    sequences.sort(key=lambda s: len(s.sequence), reverse=True)

    seen_substrings = set()
    if sequences:
        min_length = len(sequences[-1].sequence)
    
    for current_seq in sequences:
        kmers = generate_kmers(current_seq.sequence, min_length)
        if any(substring in seen_substrings for substring in generate_kmers(current_seq.sequence, min_length)):
            continue

        unique_seqs[hash_sequence(current_seq.sequence)] = current_seq
        seen_substrings.update(kmers)

    logging.info(f'Deduplicated chunk to {len(unique_seqs)} unique sequences')

    return list(unique_seqs.values())
        
def recursive_deduplication(input_file: str, file_type: str, num_threads: int) -> List[Seq]:
    # Read sequences from the temporary input file
    sequences = read_sequences_from_file(input_file, file_type)
    deduped_dict = dict()

    while True:
        if not sequences:
            break

        logging.info(f'Initial number of sequences: {len(sequences)}')
        
        total_sequences = len(sequences)
        chunk_size = max(1, total_sequences // num_threads)
        remainder = total_sequences % num_threads
        seq_chunks = []
        start_idx = 0
        for i in range(num_threads):
            end_idx = start_idx + chunk_size
            if i < remainder:
                end_idx += 1  # add one more sequence to account for remainder
            seq_chunks.append(sequences[start_idx:end_idx])
            start_idx = end_idx

        deduped_seqs = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            func = partial(deduplicate_chunk, unique_seqs=deduped_dict)
            futures = [executor.submit(func, chunk) for chunk in seq_chunks if chunk]
            concurrent.futures.wait(futures)

            for future in concurrent.futures.as_completed(futures):
                deduped_seqs.extend(future.result())

        # Add deduplicated sequences to the main set
        for seq in deduped_seqs:
            deduped_dict[hash_sequence(seq.sequence)] = seq

        # If no sequences were deduplicated in this iteration, we're done
        if len(deduped_dict) == len(sequences):
            break

        # Otherwise, shuffle the partially deduplicated sequences and repeat the process
        deduped_seqs = list(deduped_dict.values())
        random.shuffle(deduped_seqs)
        
        sequences = deduped_seqs
    
    return list(deduped_dict.values())

def deduplicate_fasta(input_file: str, file_type: str, output_file: str, num_threads: int) -> None:
    deduped_seqs = recursive_deduplication(input_file=input_file, file_type=file_type, num_threads=num_threads)

    write_sequences_to_file(deduped_seqs, output_file)
    logging.info(f'Wrote {len(deduped_seqs)} final deduplicated sequences to {output_file}')

def main():
    '''Main function that handles command-line arguments and invokes the deduplication.'''
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')

    args = parser.parse_args()

    # Validate input and output.
    if not os.path.exists(args.input_file):
        raise ValueError(f'The input file "{args.input_file}" does not exist.')
    try:
        with open(args.output_file, 'w'):
            pass
    except Exception as e:
        raise ValueError(f'Error creating output file "{args.output_file}": {e}')

    # Check file type.
    file_type = ''
    if any(ext in args.input_file for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(ext in args.input_file for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).')

    # Validate CPUs.
    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f'Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}')
        args.num_threads = min(max(args.num_threads, 1), max_threads)
    if args.num_threads >= len(read_sequences_from_file(args.input_file, file_type)) * 0.1:
        logging.warning(f'Number of sequences too low, thread count adjusted to 1')
        args.num_threads = 1

    deduplicate_fasta(args.input_file, file_type, args.output_file, args.num_threads)

if __name__ == '__main__':
    main()
