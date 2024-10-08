from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product
from functools import partial
from collections import namedtuple
from math import sqrt
import csv
import gzip
import logging
import os
from pathlib import Path
import argparse
import sys

FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']

logging.basicConfig(level=logging.INFO)

@dataclass(frozen=False)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __post_init__(self):
        object.__setattr__(self, 'sequence', self.sequence.upper())

    def __hash__(self):
        return hash((self.id, self.sequence))

    def __eq__(self, other):
        if isinstance(other, Seq):
            return self.id == other.id and self.sequence == other.sequence
        return False

    def length(self) -> int:
        return len(self.sequence)

def read_sequences(file_path: Path) -> list[Seq]:
    sequences = []
    file_type = None
    if any(file_path.suffix in ext for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(file_path.suffix in ext for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {file_path}. Expected FASTA {FASTA_EXTENSIONS} or FASTQ {FASTQ_EXTENSIONS}')
    opener = gzip.open if file_path.suffix == '.gz' else open
    if file_type == 'FASTQ':
        with opener(file_path, 'rt') as fin:
            while True:
                header_line = fin.readline().strip()
                if not header_line:
                    break
                sequence_line = fin.readline().strip()
                _ = fin.readline().strip()  # Plus line
                quality_line = fin.readline().strip()
                sequences.append(Seq(header_line[1:], sequence_line, quality_line))
    elif file_type == 'FASTA':
        with opener(file_path, 'rt') as fin:
            seq_id, seq = None, []
            for line in fin:
                line = line.strip()
                if line.startswith('>'):
                    if seq_id:
                        sequences.append(Seq(seq_id, ''.join(seq)))
                    seq_id, seq = line[1:], []
                else:
                    seq.append(line.upper())
            if seq_id:
                sequences.append(Seq(seq_id, ''.join(seq)))
    return sequences

def sequence_to_bits(dna: str) -> list:
    nucl_to_bits = {
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
    return [nucl_to_bits[nu] for nu in dna]

def reverse_complement(dna: str) -> str:
    complement = str.maketrans('ATGCRYSWKMBDHVN', 'TACGYRSWMKVHDBN')
    return dna[::-1].translate(complement)

def kmp_search(text: str, pattern: str, tolerance=0) -> list[tuple[int, int, float]]:
    def compute_lps(pattern: str) -> list[int]:
        lps = [0] * len(pattern)
        length = 0  # Length of the previous longest prefix suffix
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    lps = compute_lps(pattern)
    matches = set()
    i = 0
    j = 0
    while i < len(text) - len(pattern) + 1:
        mismatch = 0
        while j < len(pattern) and i + j < len(text):
            if text[i + j] == pattern[j]:
                j += 1
            else:
                mismatch += 1
                if mismatch > tolerance:
                    break
                j += 1
        if mismatch <= tolerance:
            sim_score = (len(pattern) - mismatch) / len(pattern)
            matches.add((i, mismatch, sim_score))
        i += 1
        j = lps[j - 1]
    return matches

def map_short_to_long(short: Seq,
                      long: Seq,
                      sim_thres: float,
                      cov_thres: float,
                      is_nucl: bool) -> list[tuple] | None:
    short, long = (long, short) if short.length() > long.length() else (short, long)
    MatchResult = namedtuple('MatchResult', ['short_id', 'long_id', 'region', 'mismatch', 'similarity', 'coverage', 'strand'])
    max_mismatch = short.length() - round(sim_thres * short.length())
    coverage = short.length() / long.length()
    if coverage < cov_thres:
        return None
    results = set()
    def process_matches(matches, strand, sequence, other_matches=None):
        for pos, mm, score in matches:
            region = f'{pos + 1}..{pos + short.length()}'
            match = MatchResult(short.id, long.id, region, mm, score, coverage, strand)
            if other_matches:
                for o_pos, o_mm, o_score in other_matches:  # Compare with reverse complement matches if available
                    if o_score > score and o_mm < mm:
                        o_region = f'{o_pos + 1}..{o_pos + short.length()}'
                        match = MatchResult(short.id, long.id, o_region, o_mm, o_score, coverage, '-' if strand == '+' else '+')
            results.add(match)
    
    if is_nucl:
        rc_short = reverse_complement(short.sequence)
        fw_matches = kmp_search(long.sequence, short.sequence, sequence_to_bits, max_mismatch)
        rc_matches = kmp_search(long.sequence, rc_short, sequence_to_bits, max_mismatch)
        if fw_matches:
            process_matches(fw_matches, '+', short.sequence, rc_matches)
        elif rc_matches:
            process_matches(rc_matches, '-', rc_short)
    else:
        matches = kmp_search(long.sequence, short.sequence, modifier=None, tolerance=max_mismatch)
        process_matches(matches, '.', short.sequence)
    formatted_results = [
        (match.short_id, match.long_id, match.region, match.mismatch, f'{match.similarity:.2f}', f'{match.coverage:.2f}', match.strand)
        for match in results
    ]
    return formatted_results if formatted_results else None

def map_sequences(query_sequences: list[Seq],
                  reference_sequences: list[Seq],
                  similarity_threshold: float,
                  coverage_threshold: float,
                  is_nucleotide: bool,
                  output_file: str,
                  num_threads: int) -> None:
    results = []
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        func = partial(map_short_to_long,
                       sim_thres=similarity_threshold,
                       cov_thres=coverage_threshold,
                       is_nucl=is_nucleotide)
        futures = [executor.submit(func, query, reference)
                   for query, reference in product(query_sequences, reference_sequences)]
        for future in as_completed(futures):
            if future.result():
                results.extend(future.result())
    headers = ['Query_ID', 'Reference_ID', 'Region', 'Mismatch', 'Similarity', 'Coverage', 'Strand']
    col_widths = {
        'Query_ID': 20,
        'Reference_ID': 20,
        'Region': 15,
        'Mismatch': 10,
        'Similarity': 10,
        'Coverage': 10,
        'Strand': 5
    }
    if results:
        results = sorted(results, key=lambda x: x[4], reverse=True)
        stdout_table(headers, results[:10], col_widths)
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)
            csv_writer.writerows(results)
    else:
        logging.info(f'No matches found.')

def truncate_string(string: str, width: int) -> str:
    return string if len(string) <= width else string[:width - 3] + '...'

def stdout_table(headers: list[str], rows: list[list[str]], col_widths: dict):
    def format_cell(value: str, width: int) -> str:
        truncated = value if len(value) <= width else value[:width - 3] + '...'
        return truncated.ljust(width)  # Left-justify to ensure even spacing

    header_row = " | ".join([format_cell(header, col_widths[header]) for header in headers])
    print(f"| {header_row} |")
    divider_row = "-" * (sum(col_widths.values()) + 3 * (len(headers) - 1) + 2)
    print(divider_row)
    for row in rows:
        formatted_row = " | ".join([format_cell(str(value), col_widths[header]) for value, header in zip(row, headers)])
        print(f"| {formatted_row} |")

def main():
    parser = argparse.ArgumentParser(description='Sequence comparison using Knuth–Morris–Pratt algorithm with Suffix Array.')
    parser.add_argument('--query', required=True, help='Path to the query input file (FASTA/FASTQ)')
    parser.add_argument('--reference', required=True, help='Path to the reference input file (FASTA/FASTQ)')
    parser.add_argument('-s', '--similarity', type=float, default=0.8,
                        help='Similarity threshold for sequence matching (default: 0.8)')
    parser.add_argument('-c', '--coverage', type=float, default=0.0,
                        help='Coverage threshold for sequence matching (default: 0.0)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file')
    parser.add_argument('-m', '--mode', choices=['nu', 'aa'], default='nu',
                        help='Comparison mode: "nu" for DNA/RNA, "aa" for proteins (default: "nu")')
    args = parser.parse_args()
    query_seqs = read_sequences(Path(args.query))
    ref_seqs = read_sequences(Path(args.reference))
    if not query_seqs or not ref_seqs:
        raise ValueError('Error: One or both input files are empty.')
    is_nucleotide = True if args.mode == 'nu' else False
    max_threads = os.cpu_count()
    if args.threads < 1 or args.threads > max_threads:
        logging.warning(f'Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}')
        args.threads = min(max(args.threads, 1), max_threads)
    if not (0 <= args.similarity < 1):
        raise ValueError('Similarity threshold should be between 0 (inclusive) and 1 (exclusive)')
    map_sequences(query_seqs, ref_seqs, args.similarity, args.coverage,
                  is_nucleotide, args.output, args.threads)

if __name__ == '__main__':
    main()
