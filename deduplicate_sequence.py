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
import tempfile


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

def generate_kmers(sequence: str, k: int) -> List[str]:
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def hash_kmers(kmers: List[str]) -> List[int]:
    return [hash_sequence(kmer) for kmer in kmers]

def is_similar_based_on_kmers(sequence1: str, sequence2: str, k: int, threshold: float) -> bool:
    # Generate k-mers for both sequences.
    kmers1 = generate_kmers(sequence1, k)
    kmers2 = generate_kmers(sequence2, k)

    # Hash the k-mers.
    hashed_kmers1 = hash_kmers(kmers1)
    hashed_kmers2 = hash_kmers(kmers2)

    # Count the occurrences of each hashed k-mer.
    counts1 = Counter(hashed_kmers1)
    counts2 = Counter(hashed_kmers2)

    # Identify the hashed k-mers present in one sequence but missing in the other.
    not_in_list2 = list((element, count) for element, count in (counts1 - counts2).items())
    not_in_list1 = list((element, count) for element, count in (counts2 - counts1).items())

    # Calculate the total number of hashed k-mers.
    total_hashed_kmers = set(counts1.keys()) | set(counts2.keys())

    # Calculate the number of differing hashed k-mers.
    count = sum(count for _, count in not_in_list1 + not_in_list2)

    return count / len(total_hashed_kmers) >= threshold if len(total_hashed_kmers) > 0 else False

# def minhash(sequence: str, k: int) -> int:
#     '''Generate a MinHash signature for a sequence using SHA-3.'''
#     min_hash_value = float('inf')

#     for i in range(len(sequence) - k + 1):
#         kmer = sequence[i:i + k]
#         hash_value = int(hashlib.sha3_256(kmer.encode()).hexdigest(), 16)
#         min_hash_value = min(min_hash_value, hash_value)

#     return min_hash_value

# def is_similar_based_on_minhash(sequence1: str, sequence2: str, k: int, threshold: float) -> bool:
#     '''Check if two sequences are similar based on their MinHash signatures.'''
#     min_hash1 = minhash(sequence1, k)
#     min_hash2 = minhash(sequence2, k)
    
#     # Calculate the similarity score (scaled between 0 and 1)
#     similarity_score = 1 - abs(min_hash1 - min_hash2) / (2**256 - 1)

#     return similarity_score

def deduplicate_chunk(sequences: List[Seq], k: int, similarity_threshold: float) -> List[Seq]:
    logging.info(f'Processing a chunk with {len(sequences)} sequences')
    sequences.sort(key=lambda s: len(s.sequence), reverse=True)  # sort by length
    unique_seqs = dict()

    for current_seq in sequences:
        # Check if a sequence is already seen
        if any(re.search(re.escape(current_seq.sequence), unique_seq.sequence) for unique_seq in unique_seqs.values()):
            continue

        # Check for similarity only if the threshold is below 1.0
        if similarity_threshold < 1.0:
            # Use non-exact match comparison
            if any(is_similar_based_on_kmers(current_seq.sequence, unique_seq.sequence, k, similarity_threshold) for unique_seq in unique_seqs.values()):
                continue

        unique_seqs[hash_sequence(current_seq.sequence)] = current_seq

    logging.info(f'Deduplicated chunk to {len(unique_seqs)} unique sequences')

    return list(unique_seqs.values())

def recursive_deduplication(input_file: str, file_type: str, num_threads: int, k: int, similarity_threshold: float) -> List[Seq]:
    deduped_dict = dict()
    temp_file_names = []

    while True:
        # Read sequences from the temporary input file
        sequences = read_sequences_from_file(input_file, file_type)

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

        temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w+t', dir='.')
        temp_file.close()

        deduped_seqs = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            func = partial(deduplicate_chunk, k=k, similarity_threshold=similarity_threshold)
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

        # Otherwise, shuffle the deduplicated sequences and write them to the input temp file
        deduped_seqs = list(deduped_dict.values())
        random.shuffle(deduped_seqs)
        
        logging.info(f'Temporarily writing {len(deduped_seqs)} sequences to {temp_file.name}')
        write_sequences_to_file(deduped_seqs, temp_file.name)

        input_file = temp_file.name
        temp_file_names.append(temp_file.name)

    # Remove all temporary files
    for temp_file in temp_file_names:
        cleanup_temp_file(temp_file)
    
    return deduped_seqs

def cleanup_temp_file(file_path: str) -> None:
    '''Remove the specified file and log if there's an error.'''
    try:
        os.remove(file_path)
    except Exception as e:
        logging.warning(f'Unable to delete temporary file {file_path}. Error: {e}')

def deduplicate_fasta(input_file: str, file_type: str, output_file: str, num_threads: int, k: int, similarity_threshold: float) -> None:
    deduped_seqs = recursive_deduplication(input_file=input_file, file_type=file_type, num_threads=num_threads, k=k, similarity_threshold=similarity_threshold)
    
    # Final scan using deduplicate_chunk to ensure all sequences are deduplicated
    logging.info(f'Final scan to deduplicate {len(deduped_seqs)} sequences')
    deduped_seqs_final = deduplicate_chunk(deduped_seqs, k, similarity_threshold)

    write_sequences_to_file(deduped_seqs_final, output_file)
    logging.info(f'Wrote {len(deduped_seqs_final)} final deduplicated sequences to {output_file}')

def main():
    '''Main function that handles command-line arguments and invokes the deduplication.'''
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')
    parser.add_argument('-k', '--k_mer', type=int, default=59, help='Number of k to generate k-mers. Default is 59.')
    parser.add_argument('-s', '--similarity_threshold', type=float, default=1.0, help='Similarity threshold based on LCS. Default is 1.0 (exact matches).')

    args = parser.parse_args()

    # Validate input and output.
    if not os.path.exists(args.input_file):
        raise ValueError(f"The input file '{args.input_file}' does not exist.")
    try:
        with open(args.output_file, 'w'):
            pass
    except Exception as e:
        raise ValueError(f"Error creating output file '{args.output_file}': {e}")

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

    deduplicate_fasta(args.input_file, file_type, args.output_file, args.num_threads, args.k_mer, args.similarity_threshold)

if __name__ == '__main__':
    main()
