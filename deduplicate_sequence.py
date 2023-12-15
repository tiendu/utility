import concurrent.futures
import argparse
import os
import gzip
import logging
import hashlib
from itertools import groupby
from dataclasses import dataclass
from typing import List, Set, Tuple
import tempfile

# Setup logging to display messages with INFO level and above.
logging.basicConfig(level=logging.INFO)

# Data class to represent a sequence.
@dataclass
class Seq:
    id: str
    sequence: str

def read_sequences_from_file(file_path: str) -> List[Seq]:
    sequences = []
    
    # Determine the appropriate file opening mode based on the file extension.
    opener = gzip.open if file_path.endswith('.gz') else open
    
    with opener(file_path, 'rt') as fin:
        faiter = (x[1] for x in groupby(fin, lambda line: line[0] == ">"))
        for header in faiter:
            headerstr = header.__next__().strip().replace('>', '')
            name = headerstr.split()[0]
            seq = "".join(s.strip().upper() for s in faiter.__next__())
            sequences.append(Seq(name, seq))
    return sequences

def hash_sequence(sequence: str) -> str:
    """Returns a SHA-256 hash of a sequence."""
    return hashlib.sha256(sequence.encode()).hexdigest()

def deduplicate_chunk(sequences: List[Seq]) -> Tuple[List[Seq], Set[str]]:
    """Deduplicates a chunk of sequences."""
    logging.info(f"Processing a chunk with {len(sequences)} sequences.")
    sequences.sort(key=lambda s: len(s.sequence), reverse=True)  # sort by length
    unique_seqs = []
    local_seen_hashes = set()

    for current_seq in sequences:
        seq_hash = hash_sequence(current_seq.sequence)

        # Check for exact match.
        if seq_hash in local_seen_hashes:
            continue

        # Check if the sequence is contained within any of the unique sequences.
        is_contained = any(current_seq.sequence in unique_seq.sequence for unique_seq in unique_seqs)
        if not is_contained:
            unique_seqs.append(current_seq)
            local_seen_hashes.add(seq_hash)
    logging.info(f"Deduplicated chunk to {len(unique_seqs)} unique sequences.")
    return unique_seqs, local_seen_hashes

def parallel_deduplication(sequences: List[Seq], num_threads: int) -> List[Seq]:
    """Deduplicates sequences using multiple threads."""
    logging.info(f"Beginning parallel deduplication with {num_threads} threads.")
    chunk_size = len(sequences) // num_threads
    seq_chunks = [sequences[i:i + chunk_size] for i in range(0, len(sequences), chunk_size)]

    all_deduped_chunks = []
    global_seen_hashes = set()

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(deduplicate_chunk, seq_chunk) for seq_chunk in seq_chunks]
        
        for future in concurrent.futures.as_completed(futures):
            deduped_chunk, chunk_seen_hashes = future.result()
            all_deduped_chunks.extend(deduped_chunk)
            global_seen_hashes.update(chunk_seen_hashes)

    # Deduplicate the consolidated results from all chunks
    final_results, _ = deduplicate_chunk(all_deduped_chunks)
    logging.info(f"Consolidated all chunks to {len(final_results)} unique sequences.")
    return final_results

def deduplicate_fasta(input_file: str, output_file: str, num_threads: int) -> None:
    """Main deduplication function. Reads sequences from the input file, deduplicates them, and writes the unique sequences to the output file."""
    sequences = read_sequences_from_file(input_file)
    result = parallel_deduplication(sequences, num_threads)

    with open(output_file, 'w') as outfile:
        for seq_obj in result:
            outfile.write(">" + seq_obj.id + "\n")
            outfile.write(seq_obj.sequence + "\n")
    logging.info(f"Wrote {len(result)} unique sequences to {output_file}.")

def main():
    """Main function that handles command-line arguments and invokes the deduplication."""
    parser = argparse.ArgumentParser(description='Deduplicate FASTA sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input FASTA file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output FASTA file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')

    args = parser.parse_args()

    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Adjusting thread count to be between 1 and {max_threads}.")
        args.num_threads = min(max(args.num_threads, 1), max_threads)

    deduplicate_fasta(args.input_file, args.output_file, args.num_threads)

if __name__ == '__main__':
    main()
