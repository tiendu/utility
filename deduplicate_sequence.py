import argparse
import os
import gzip
import logging
import hashlib
import random
from functools import partial
from itertools import groupby
from dataclasses import dataclass
from typing import List, Generator, Dict, Set, Callable, Any
import concurrent.futures

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna']
LINE_LENGTH = 80
CHUNK_SIZE = 100_000

# Setup logging to display messages with INFO level and above.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''

    def __hash__(self) -> int:
        return hash((self.id, self.sequence, self.quality))

    def length(self) -> int:
        return len(self.sequence)

def read_sequences_from_file(file_path: str, file_type: str) -> List[Seq]:
    sequences = []
    opener = gzip.open if file_path.endswith('.gz') else open

    if file_type == 'FASTQ':
        with opener(file_path, 'rt') as fin:
            groups = groupby(enumerate(fin), key=lambda x: x[0] // 4)
            for _, group in groups:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                name = header_line[1:]
                seq = sequence_line.upper()
                qual = quality_line.upper()
                sequences.append(Seq(name, seq, qual))
    elif file_type == 'FASTA':
        with opener(file_path, 'rt') as fin:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == '>'))
            for header in faiter:
                headerstr = next(header).strip()
                name = headerstr[1:]
                seq = ''.join(s.strip().upper() for s in next(faiter))
                sequences.append(Seq(name, seq))
    return sequences

def write_sequences_to_file(sequences: List[Seq], file_path: str) -> None:
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'wt') as f:
        for seq in sequences:
            if seq.quality == '':
                f.write(f'>{seq.id}\n')
                for i in range(0, len(seq.sequence), LINE_LENGTH):
                    f.write(seq.sequence[i:i+LINE_LENGTH] + '\n')
            else:
                f.write(f'@{seq.id}\n{seq.sequence}\n+\n{seq.quality}\n')

def generate_kmers(string: str, k: int) -> Generator[str, None, None]:
    for i in range(len(string) - k + 1):
        yield string[i:i+k]

def hash_string(string: str, hash_function=hashlib.sha3_256) -> str:
    return hash_function(string.encode()).hexdigest()

def round_robin_divide(items: List[Any], chunk_size: int, num_threads: int, key: Callable[[Any], Any], is_descending: bool) -> List[List[Any]]:
    def is_sorted(lst: List[Any], comparison_func: Callable[[Any, Any], bool]) -> bool:
        return all(comparison_func(lst[i], lst[i + 1]) for i in range(len(lst) - 1))

    if is_descending:
        order = lambda x, y: x >= y
    else:
        order = lambda x, y: x <= y

    if not is_sorted([key(item) for item in items], order):
        items = sorted(items, key=key, reverse=is_descending)

    chunks = [[] for _ in range(num_threads)]
    for i, item in enumerate(items):
        chunks[i % num_threads].append(item)

    divided_chunks = []
    for chunk in chunks:
        for i in range(0, len(chunk), chunk_size):
            divided_chunks.append(chunk[i:i + chunk_size])

    return divided_chunks

def deduplicate_chunk(sequences: List[Seq], uniq_kmers: Set[str], min_length: int) -> List[Seq]:
    local_uniq_seqs = {}
    local_uniq_kmers = set()

    for sequence in sequences:
        kmer_hashes = set(hash_string(kmer) for kmer in generate_kmers(sequence.sequence, min_length))
        if all(kmer_hash in uniq_kmers or kmer_hash in local_uniq_kmers for kmer_hash in kmer_hashes):
            continue
        local_uniq_seqs[hash_string(sequence.sequence)] = sequence
        local_uniq_kmers.update(kmer_hashes)

    uniq_kmers.update(local_uniq_kmers)

    return list(local_uniq_seqs.values())

def deduplicate_concurrently(sequences: List[Seq], num_threads: int) -> List[Seq]:
    chunk_size = max(1, min(CHUNK_SIZE, len(sequences) // num_threads))

    while sequences:
        logging.info(f'Current number of sequences: {len(sequences)}')
        shared_sequences: Set[Seq] = set()
        shared_kmers: Set[str] = set()
        sequences = sorted(sequences, key=lambda sequence: sequence.length(), reverse=True)
        min_length = sequences[-1].length()
        chunks = round_robin_divide(sequences, chunk_size, num_threads, key=lambda sequence: sequence.length(), is_descending=True)

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
            func = partial(deduplicate_chunk, uniq_kmers=shared_kmers, min_length=min_length)
            futures = [executor.submit(func, chunk) for chunk in chunks]
            concurrent.futures.wait(futures)
            for future in concurrent.futures.as_completed(futures):
                shared_sequences.update(future.result())

        if len(shared_sequences) == len(sequences):
            break
        sequences = list(shared_sequences)
        random.shuffle(sequences)

    logging.info(f'Final number of sequences: {len(shared_sequences)}')
    return list(shared_sequences)

def main() -> None:
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        raise ValueError(f'The input file "{args.input_file}" does not exist.')
    try:
        with open(args.output_file, 'w'):
            pass
    except Exception as e:
        raise ValueError(f'Error creating output file "{args.output_file}": {e}')

    file_type = ''
    if any(ext in args.input_file for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(ext in args.input_file for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).')

    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f'Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}')
        args.num_threads = min(max(args.num_threads, 1), max_threads)

    sequences = read_sequences_from_file(args.input_file, file_type)
    if not sequences:
        raise ValueError('No sequences detected!')

    if args.num_threads >= len(sequences) * 0.1:
        logging.warning(f'Number of sequences too low, thread count adjusted to 1')
        args.num_threads = 1

    deduplicated_sequences = deduplicate_concurrently(sequences, args.num_threads)
    write_sequences_to_file(deduplicated_sequences, args.output_file)

if __name__ == '__main__':
    main()
