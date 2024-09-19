from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product
from functools import partial
from math import sqrt
import csv
import gzip
from pathlib import Path
import argparse
import sys

FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']

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

def kmerize(string: str, k: int, modifier=None) -> list[str]:
    if len(string) < k:
        return []
    if modifier:
        return [modifier(string[i:i+k]) for i in range(len(string) - k + 1)]
    else:
        return [string[i:i+k] for i in range(len(string) - k + 1)]

def compare_kmer_sets(kmers1: set, kmers2: set) -> int:
    def compare_kmers(kmer_bits1: set, kmer_bits2: set):
        return all(bit1 & bit2 for bit1, bit2 in zip(kmer_bits1, kmer_bits2))

    if not kmers1 or not kmers2:
        return 0
    total_matches = 0
    for kmer_bits1 in kmers1:
        best_match = max(compare_kmers(kmer_bits1, kmer_bits2) for kmer_bits2 in kmers2)
        total_matches += best_match
    return total_matches

def calculate_kmer_similarity(seq1: str, seq2: str, k: int) -> float:
    kmers1_bits = kmerize(seq1, k, sequence_to_bits)
    kmers2_bits = kmerize(seq2, k, sequence_to_bits)
    total_matches = compare_kmer_sets(kmers1_bits, kmers2_bits)
    total_kmers = max(len(kmers1_bits), len(kmers2_bits))
    similarity = total_matches / total_kmers if total_kmers > 0 else 0
    return similarity

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

def short_to_sub_long(i: int,
                      k: int,
                      short: Seq,
                      long: Seq,
                      sim_thresh: float,
                      cov_thresh: float,
                      is_nucleotide: bool,
                      is_circular: bool) -> tuple[str, str, str, float, float, str] | None:
    def get_similarity_and_strand(fw_match: str, rc_match: str = '') -> tuple[float, str]:
        if is_nucleotide:
            rc_sim = calculate_kmer_similarity(short.sequence, rc_match, k)
            fw_sim = calculate_kmer_similarity(short.sequence, fw_match, k)
            if fw_sim >= rc_sim:
                return fw_sim, '+'
            return rc_sim, '-'
        return calculate_kmer_similarity(short.sequence, fw_match, k), '.'

    fw_match = long.sequence[i:i + len(short.sequence)]
    rc_match = reverse_complement(fw_match) if is_nucleotide else ''
    best_sim, strand = get_similarity_and_strand(fw_match, rc_match)
    if is_circular:
        circular_long = long.sequence + long.sequence
        fw_match_circ = circular_long[i:i + len(short.sequence)]
        rc_match_circ = reverse_complement(fw_match_circ) if is_nucleotide else ''
        best_sim_circ, strand_circ = get_similarity_and_strand(fw_match_circ, rc_match_circ)
        if best_sim_circ > best_sim:
            best_sim = best_sim_circ
            strand = strand_circ
            match_type = 'circular'
        else:
            match_type = 'linear'
    else:
        match_type = 'linear'
    coverage = len(fw_match) / len(long.sequence)
    if coverage >= cov_thresh and best_sim >= sim_thresh:
        start = i + 1
        end = (i + len(short.sequence)) % len(long.sequence) if is_circular else i + len(short.sequence)
        position = f'{start}..{end}' if match_type == 'linear' else f'{start}..{end} (circular)'
        return short.id, long.id, position, round(best_sim, 3), round(coverage, 3), strand
    return None

def map_short_to_long(short: Seq,
                      long: Seq,
                      similarity_threshold: float,
                      coverage_threshold: float,
                      is_nucleotide: bool,
                      is_circular: bool=False,
                      num_threads: int=1) -> list[tuple[str, str, str, float, float, str]]:
    if len(short.sequence) > len(long.sequence):
        short, long = long, short 
    k = max(len(short.sequence) // 5 | 1, 3)
    results = []
    indices = range(len(long.sequence) - len(short.sequence) + 1)
    if is_circular:
        indices = range((len(long.sequence) - len(short.sequence) + 1) * 2)
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        func = partial(short_to_sub_long, k=k,
                       short=short, long=long,
                       sim_thresh=similarity_threshold, cov_thresh=coverage_threshold,
                       is_nucleotide=is_nucleotide, is_circular=is_circular)
        
        futures = [executor.submit(func, i) for i in indices]
        
        for future in as_completed(futures):
            result = future.result()
            if result:
                results.append(result)
                
    return results

def map_sequences(query_sequences: list[Seq],
                  reference_sequences: list[Seq],
                  similarity_threshold: float,
                  coverage_threshold: float,
                  is_nucleotide: bool,
                  is_circular: bool,
                  output_file: str,
                  num_threads: int=1) -> None:
    results = []
    for query, reference in product(query_sequences, reference_sequences):
        result = map_short_to_long(query, reference,
                                   similarity_threshold, coverage_threshold,
                                   is_nucleotide, is_circular,
                                   num_threads)
        results.extend(result)
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

def main():
    parser = argparse.ArgumentParser(description='Sequence comparison using k-mers.')
    parser.add_argument('--query', required=True, help='Path to the query input file (FASTA/FASTQ)')
    parser.add_argument('--reference', required=True, help='Path to the reference input file (FASTA/FASTQ)')
    parser.add_argument('-s', '--similarity', type=float, default=0.8,
                        help='Similarity threshold for sequence matching (default: 0.8)')
    parser.add_argument('-c', '--coverage', type=float, default=0.0,
                        help='Coverage threshold for sequence matching (default: 0.0)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file')
    parser.add_argument('--mode', choices=['nu', 'aa'], default='nu',
                        help="Comparison mode: 'nu' for DNA/RNA, 'aa' for proteins (default: 'nu')")
    parser.add_argument('--topology', choices=['linear', 'circular'], default='linear',
                        help='Determine whether a sequence is linear/circular')
    args = parser.parse_args()
    query_seqs = read_sequences(Path(args.query))
    ref_seqs = read_sequences(Path(args.reference))
    if not query_seqs or not ref_seqs:
        print('Error: One or both input files are empty.')
        sys.exit(1)
    is_circular = False if args.topology == 'linear' else True
    is_nucleotide = True if args.mode == 'nu' else False
    if not is_nucleotide and is_circular:
        print("Error: Amino acid sequences don't have circular topology")
        sys.exit(1)
    map_sequences(query_seqs, ref_seqs, args.similarity, args.coverage,
                  is_nucleotide, is_circular, args.output, args.threads)

if __name__ == '__main__':
    main()
