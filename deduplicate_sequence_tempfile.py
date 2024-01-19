import argparse
import os
import gzip
import logging
import hashlib
import random
from collections import Counter
from multiprocessing import Manager, Lock
from functools import partial, lru_cache
from itertools import groupby, product
from dataclasses import dataclass
import concurrent.futures
from typing import List, Set, Tuple, Dict
import tempfile


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

    if file_type == "FASTQ":
        with opener(file_path, 'rt') as fin:
            groups = groupby(enumerate(fin), key=lambda x: x[0] // 4)
            for _, group in groups:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                name = header_line[1:]  # Removing "@" character
                seq = sequence_line.upper()
                qual = quality_line.upper()
                sequences.append(Seq(name, seq, qual))
    elif file_type == "FASTA":
        with opener(file_path, 'rt') as fin:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == ">"))
            for header in faiter:
                headerstr = next(header).strip()
                name = headerstr[1:]  # Removing ">" character
                seq = "".join(s.strip().upper() for s in next(faiter))
                sequences.append(Seq(name, seq))
    return sequences

def write_sequences_to_file(sequences: List[Seq], file_path: str) -> None:
    with open(file_path, 'w') as f:
        for seq in sequences:
            if seq.quality == '':
                f.write(f">{seq.id}\n")
                # Write sequence with a maximum line length of 80 characters
                for i in range(0, len(seq.sequence), 80):
                    f.write(seq.sequence[i:i+80] + "\n")
            else:
                f.write(f"@{seq.id}\n{seq.sequence}\n+\n{seq.quality}\n")

@lru_cache(maxsize=None)
def hash_sequence(sequence: str, hash_function=hashlib.sha3_256) -> str:
    """Return the hash of a sequence using the specified hash function."""
    return hash_function(sequence.encode()).hexdigest()

@lru_cache(maxsize=None)
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

@lru_cache(maxsize=None)
def minhash(sequence: str, k: int) -> int:
    """Generate a MinHash signature for a sequence using SHA-3."""
    min_hash_value = float('inf')

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        hash_value = int(hashlib.sha3_256(kmer.encode()).hexdigest(), 16)
        min_hash_value = min(min_hash_value, hash_value)

    return min_hash_value

def is_similar_based_on_minhash(sequence1: str, sequence2: str, k: int, threshold: float) -> bool:
    """Check if two sequences are similar based on their MinHash signatures."""
    min_hash1 = minhash(sequence1, k)
    min_hash2 = minhash(sequence2, k)
    
    # Calculate the similarity score (scaled between 0 and 1)
    similarity_score = 1 - abs(min_hash1 - min_hash2) / (2**256 - 1)
    
    return similarity_score

def deduplicate_chunk(sequences: List[Seq], k: int, similarity_threshold: float, seen_hashes: Dict[str, bool], lock: Lock) -> List[Seq]:
    logging.info(f"Processing a chunk with {len(sequences)} sequences.")
    sequences.sort(key=lambda s: len(s.sequence), reverse=True)  # sort by length
    unique_seqs = []

    for current_seq in sequences:
        seq_hash = hash_sequence(current_seq.sequence)

        with lock:
            if seen_hashes.get(seq_hash, False):
                continue

            seen_hashes[seq_hash] = True

            if similarity_threshold == 1.0:
                # Check if a sequence is already seen using its hash
                if any(current_seq.sequence in unique_seq.sequence for unique_seq in unique_seqs):
                    continue
            else:
                # Check if a sequence is contained in any of the unique sequences
                if any(is_similar_based_on_kmers(current_seq.sequence, unique_seq.sequence, k, similarity_threshold) for unique_seq in unique_seqs):
                    continue

        unique_seqs.append(current_seq)

    logging.info(f"Deduplicated chunk to {len(unique_seqs)} unique sequences.")
    return unique_seqs

def recursive_deduplication(input_file: str, file_type: str, num_threads: int, k: int, similarity_threshold: float, seen_hashes_shared: Manager().dict, lock: Lock) -> List[Seq]:
    deduped_seqs_all = set()
    temp_file_names = []

    while True:
        # Read sequences from the temporary input file
        sequences = read_sequences_from_file(input_file, file_type)

        if not sequences:
            break

        logging.info(f"Initial number of sequences: {len(sequences)}")
        
        total_sequences = len(sequences)
        chunk_size = max(1, total_sequences // num_threads)
        seq_chunks = [sequences[i:i + chunk_size] for i in range(0, total_sequences, chunk_size)]

        temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w+t', dir='.')
        temp_file.close()

        deduped_seqs = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            func = partial(deduplicate_chunk, k=k, similarity_threshold=similarity_threshold, seen_hashes=seen_hashes_shared, lock=lock)
            futures = [executor.submit(func, chunk) for chunk in seq_chunks if chunk]

            for future in concurrent.futures.as_completed(futures):
                deduped_seqs.extend(future.result())

        # Add deduplicated sequences to the main set
        deduped_seqs_all.update(set(deduped_seqs))

        # If no sequences were deduplicated in this iteration, we're done
        if len(deduped_seqs_all) == len(sequences):
            break

        # Otherwise, shuffle the deduplicated sequences and write them to the input temp file
        sequences = list(deduped_seqs_all)
        random.shuffle(sequences)
        
        logging.info(f"Temporarily writing {len(sequences)} sequences to {temp_file.name}.")
        write_sequences_to_file(sequences, temp_file.name)

        input_file = temp_file.name
        temp_file_names.append(temp_file.name)

    # Remove all temporary files
    for temp_file in temp_file_names:
        cleanup_temp_file(temp_file)
    
    return list(deduped_seqs_all)

def cleanup_temp_file(file_path: str) -> None:
    """Remove the specified file and log if there's an error."""
    try:
        os.remove(file_path)
    except Exception as e:
        logging.warning(f"Unable to delete temporary file {file_path}. Error: {e}")

def deduplicate_fasta(input_file: str, file_type: str, output_file: str, num_threads: int, k: int, similarity_threshold: float) -> None:
    with Manager() as manager:
        seen_hashes_shared = manager.dict()
        lock = manager.Lock()
        deduped_seqs = recursive_deduplication(input_file=input_file, file_type=file_type, num_threads=num_threads, k=k, similarity_threshold=similarity_threshold, seen_hashes_shared=seen_hashes_shared, lock=lock)
    
    write_sequences_to_file(deduped_seqs, output_file)
    logging.info(f"Wrote {len(deduped_seqs)} final deduplicated sequences to {output_file}")

def main():
    """Main function that handles command-line arguments and invokes the deduplication."""
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')
    parser.add_argument('-k', '--k_mer', type=int, default=59, help='Number of k to generate k-mers. Default is 59.')
    parser.add_argument('-s', '--similarity_threshold', type=float, default=1.0, help='Similarity threshold based on LCS. Default is 1.0 (exact matches).')

    args = parser.parse_args()

    file_type = ''
    if any(ext in args.input_file for ext in ['.fastq', '.fq']):
        file_type = "FASTQ"
    elif any(ext in args.input_file for ext in ['.fasta', '.fa', '.fna']):
        file_type = "FASTA"
    else:
        raise ValueError(f"Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).")

    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}.")
        args.num_threads = min(max(args.num_threads, 1), max_threads)
    if args.num_threads >= len(read_sequences_from_file(args.input_file, file_type)) * 0.1:
        logging.warning(f"Number of sequences too low. Adjusting thread count to 1.")
        args.num_threads = 1

    deduplicate_fasta(args.input_file, file_type, args.output_file, args.num_threads, args.k_mer, args.similarity_threshold)

if __name__ == '__main__':
    main()
