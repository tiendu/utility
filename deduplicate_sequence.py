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
from concurrent.futures import ThreadPoolExecutor, wait, as_completed
from collections import defaultdict

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']
LINE_LENGTH = 80

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

def hash_string(string: str, hash_function=hashlib.sha3_256) -> str:
    return hash_function(string.encode()).hexdigest()

def kmerize(string: str, k: int, modifier=None) -> list[str]:
    if len(string) < k:
        return []
    if modifier:
        return [modifier(string[i:i+k]) for i in range(len(string) - k + 1)]
    return [string[i:i+k] for i in range(len(string) - k + 1)]

def compare_two_dnas(dna1: str, dna2: str) -> bool:
    def dna_to_bits(dna: str) -> list:
        nu_to_bits = {
            'A': 0b0001,  # A = 1
            'T': 0b0010,  # T = 2
            'G': 0b0100,  # G = 4
            'C': 0b1000,  # C = 8
            'W': 0b0011,  # W = A | T = 3
            'R': 0b0101,  # R = A | G = 5
            'Y': 0b1010,  # Y = C | T = 10
            'S': 0b1100,  # S = G | C = 12
            'K': 0b0110,  # K = G | T = 6
            'M': 0b1001,  # M = A | C = 9
            'B': 0b1110,  # B = C | G | T = 14
            'D': 0b0111,  # D = A | G | T = 7
            'H': 0b1011,  # H = A | C | T = 11
            'V': 0b1101,  # V = A | C | G = 13
            'N': 0b1111   # N = A | C | G | T = 15 (any nucleotide)
        }
        return [nu_to_bits[nu.upper()] for nu in dna]
    
    def compare_bits(bits1: list, bits2: list) -> bool:
        return all(bit1 & bit2 for bit1, bit2 in zip(bits1, bits2))

    long = dna1 if len(dna1) > len(dna2) else dna2
    short = dna1 if len(dna1) <= len(dna2) else dna2
    long_bits = dna_to_bits(long)
    short_bits = dna_to_bits(short)
    for i in range(len(long) - len(short) + 1):
        if compare_bits(long_bits[i:i + len(short)], short_bits):
            return True
    return False

def deduplicate_chunk(seqs: List[Seq], mode: str, k: int=5) -> Dict[str, Seq]:
    seqs = sorted(seqs, key=lambda seq: seq.length(), reverse=True)
    uniq_seqs = {}
    kmer_dict = defaultdict(list)
    def compare_with_kmer_seqs(seq: Seq, compare_func, k: int) -> bool:
        seq_len = len(seq.sequence)
        # Case 1: Sequence is shorter than k
        if seq_len < k:
            # Check if the short sequence is a subsequence of any k-mer in the dictionary
            for kmer, seq_list in kmer_dict.items():
                if compare_func(seq.sequence, kmer):  # Compare short sequence with the k-mer
                    # Compare the short sequence with all sequences associated with the k-mer
                    for existing_seq in seq_list:
                        if compare_func(existing_seq.sequence, seq.sequence):
                            return True
            return False
        # Case 2: Sequence length is >= k
        seq_kmers = kmerize(seq.sequence, k)  # Generate k-mers from the sequence
        for kmer in seq_kmers:
            for dict_kmer in kmer_dict:  # Check if any k-mer from the dictionary matches
                if compare_func(kmer, dict_kmer):
                    for existing_seq in kmer_dict[dict_kmer]:
                        if compare_func(existing_seq.sequence, seq.sequence):
                            return True
        return False

    if mode == 'nu':  # Nucleotide comparison
        compare_func = compare_two_dnas
    elif mode == 'aa':  # Amino acid comparison
        compare_func = lambda s1, s2: s1 in s2 or s2 in s1
    for seq in seqs:
        if not compare_with_kmer_seqs(seq, compare_func, k):
            uniq_seqs[hash_string(seq.sequence)] = seq
            seq_kmers = kmerize(seq.sequence, k)
            for kmer in seq_kmers:
                kmer_dict[kmer].append(seq)
    return uniq_seqs

def deduplicate_concurrently(sequences: List[Seq], num_threads: int, mode: str) -> List[Seq]:
    deduped_seqs = dict()
    while sequences:
        logging.info(f'Current number of sequences: {len(sequences)}')
        size = max(1, len(sequences) // num_threads)  # Calculate chunk size
        chunks = [sequences[i:i + size] for i in range(0, len(sequences), size)]  # Proper chunking
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            func = partial(deduplicate_chunk,
                           mode=mode)
            futures = [executor.submit(func, chunk) for chunk in chunks]
            wait(futures)
            for future in as_completed(futures):
                deduped_seqs.update(future.result())
        if len(deduped_seqs) == len(sequences):
            break
        sequences = list(deduped_seqs.values())
        random.shuffle(sequences)
    sequences = list(deduplicate_chunk(sequences, mode).values())
    logging.info(f'Final number of sequences: {len(sequences)}')
    return sequences

def main() -> None:
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')
    parser.add_argument('-m', '--mode', choices=['nu', 'aa'], default='nu',
                        help="Comparison mode: 'nu' for DNA/RNA, 'aa' for proteins (default: 'nu')")
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
        raise ValueError(f'Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna, .faa) or FASTQ (.fastq, .fq).')
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
    deduplicated_sequences = deduplicate_concurrently(sequences, args.num_threads, args.mode)
    write_sequences_to_file(deduplicated_sequences, args.output_file)

if __name__ == '__main__':
    main()
