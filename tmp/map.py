from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
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
    return [nucl_to_bits[nu.upper()] for nu in dna]

def reverse_complement(dna: str) -> str:
    complement = str.maketrans('ATGCRYSWKMBDHVN', 'TACGYRSWMKVHDBN')
    return dna[::-1].translate(complement)

def build_suffix_array(string: str) -> list:
    suffixes = [(string[i:], i) for i in range(len(string))]
    suffixes = sorted(suffixes, key=lambda substring: substring[0], reverse=False)
    return suffixes

def compute_lps(pattern: list) -> list:
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

def kmp_search(text: str, pattern: str, modifier=None, tolerance=0) -> list[tuple[int, float]]:
    def match(a, b):
        if isinstance(a, int) and isinstance(b, int):
            return (a & b) != 0
        return a == b

    if modifier:
        text = modifier(text)
        pattern = modifier(pattern)
    else:
        text = list(text)
        pattern = list(pattern)
    if tolerance >= len(pattern):
        return set()
    lps = compute_lps(pattern)
    occurrences = set()
    for i in range(len(text)):
        mismatch = 0
        j = 0
        while j < len(pattern) and (i + j) < len(text):
            if match(text[i + j], pattern[j]):
                j += 1
            else:
                mismatch += 1
                if mismatch > tolerance:
                    break
                j += 1
        if j == len(pattern) and mismatch <= tolerance:
            sim_score = (len(pattern) - mismatch) / len(pattern)
            occurrences.add((i, sim_score))
    return occurrences

def search_with_suffix_and_kmp(text: str, pattern: str, modifier=None, threshold: float=0.0) -> set:
    results = set()
    suffixes = build_suffix_array(text)
    tolerance = len(pattern) - round(threshold * len(pattern))
    for suffix, index in suffixes:
        occurences = kmp_search(suffix, pattern, modifier, tolerance)
        if occurences:
            for occurence in occurences:
                adj_pos = index + occurence[0]
                sim_score = occurence[1]
                results.add((adj_pos, sim_score))
    return results

def map_short_to_long(short: Seq,
                      long: Seq,
                      sim_thres: float,
                      cov_thres: float,
                      is_nucl: bool) -> list[tuple[str, str, str, float, float, str]] | None:
    short, long = (long, short) if short.length() > long.length() else (short, long)
    MatchResult = namedtuple('MatchResult', ['short_id', 'long_id', 'region', 'similarity', 'coverage', 'strand'])
    coverage = short.length() / long.length()
    if coverage < cov_thres:
        return None  # Exit early if coverage is below the threshold
    results = []
    covered_regions = set()
    if is_nucl:
        rc_short = reverse_complement(short.sequence)
        fw_matches = search_with_suffix_and_kmp(long.sequence, short.sequence, sequence_to_bits, sim_thres)
        rc_matches = search_with_suffix_and_kmp(long.sequence, rc_short, sequence_to_bits, sim_thres)
        if fw_matches:
            for fw_pos, fw_score in fw_matches:
                fw_region = f'{fw_pos + 1}..{fw_pos + len(short.sequence)}'
                best_match_fw = MatchResult(short.id, long.id, fw_region, fw_score, coverage, '+')
                results.append(best_match_fw)
                covered_regions.add((fw_pos, fw_pos + len(short.sequence)))  # Mark this region as covered by a forward match
        if rc_matches:
            for rc_pos, rc_score in rc_matches:
                if all(not (rc_pos <= fw_end and rc_pos + len(short.sequence) >= fw_start)
                       for fw_start, fw_end in covered_regions):  # Check if this reverse match region overlaps with any forward match
                    rc_region = f'{rc_pos + 1}..{rc_pos + len(short.sequence)}'
                    best_match_rc = MatchResult(short.id, long.id, rc_region, rc_score, coverage, '-')
                    results.append(best_match_rc)
    else:
        matches = search_with_suffix_and_kmp(long.sequence, short.sequence, modifier=None, threshold=sim_thres)
        for pos, score in matches:
            region = f'{pos + 1}..{pos + len(short.sequence)}'
            results.append(MatchResult(short.id, long.id, region, score, coverage, '.'))
    formatted_results = [
        (match.short_id, match.long_id, match.region, f'{match.similarity:.2f}', f'{match.coverage:.2f}', match.strand)
        for match in results
    ]
    return formatted_results if formatted_results else None

def map_sequences(query_sequences: list[Seq],
                  reference_sequences: list[Seq],
                  similarity_threshold: float,
                  coverage_threshold: float,
                  is_nucleotide: bool,
                  output_file: str,
                  num_threads: int=1) -> None:
    results = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
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
    headers = ['Query_ID', 'Reference_ID', 'Region', 'Similarity', 'Coverage', 'Strand']
    col_widths = {
        'Query_ID': 20,
        'Reference_ID': 20,
        'Region': 15,
        'Similarity': 10,
        'Coverage': 10,
        'Strand': 5
    }
    if results:
        results = sorted(results, key=lambda x: x[3], reverse=True)
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
                        help="Comparison mode: 'nu' for DNA/RNA, 'aa' for proteins (default: 'nu')")
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
