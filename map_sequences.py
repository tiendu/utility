from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product, islice, groupby
from functools import partial
import csv
import gzip
import logging
from pathlib import Path
import argparse

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

    def length(self) -> int:
        return len(self.sequence)

def read_sequences(file_path: str) -> list[Seq]:
    sequences = []
    file_type = (
        'FASTQ' if any(ext in file_path for ext in FASTQ_EXTENSIONS)
        else 'FASTA' if any(ext in file_path for ext in FASTA_EXTENSIONS)
        else ValueError(f'Unrecognized file extension for {file_path}.')
    )
    opener = gzip.open if file_path.endswith('.gz') else open
    file_handler = opener(file_path, 'rt')
    if file_type == 'FASTQ':
        groups = groupby(enumerate(file_handler), key=lambda x: x[0] // 4)
        for _, group in groups:
            try:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                seqid = header_line[1:]
                seq = sequence_line
                qual = quality_line
                sequences.append(Seq(seqid, seq, qual))
            except ValueError:
                logging.warning('Malformed FASTQ entry encountered. Skipping')
    elif file_type == 'FASTA':
        faiter = (x[1] for x in groupby(file_handler, lambda line: line[0] == '>'))
        for header in faiter:
            try:
                header_str = next(header).strip()
                seqid = header_str[1:]
                seq = ''.join(s.strip() for s in next(faiter))
                sequences.append(Seq(seqid, seq))
            except StopIteration:
                logging.warning('Incomplete FASTA entry encountered. Skipping')
    file_handler.close()
    return sequences

def dna_to_bits(dna: str) -> list:
    nucl_to_bits = {
        'A': 0b0001, 'T': 0b0010, 'G': 0b0100, 'C': 0b1000,
        'W': 0b0011, 'R': 0b0101, 'Y': 0b1010, 'S': 0b1100,
        'K': 0b0110, 'M': 0b1001, 'B': 0b1110, 'D': 0b0111,
        'H': 0b1011, 'V': 0b1101, 'N': 0b1111
    }
    return [nucl_to_bits[nu] for nu in dna]

def reverse_complement(dna: str) -> str:
    complement = str.maketrans('ATGCRYSWKMBDHVN', 'TACGYRSWMKVHDBN')
    return dna[::-1].translate(complement)

def search(text: str, pattern: str, modifier=None, tolerance=0) -> set:
    text = modifier(text) if modifier else text
    pattern = tuple(modifier(pattern)) if modifier else pattern
    matches = set()
    for i in range(len(text) - len(pattern) + 1):
        mismatch = sum(1 for j in range(len(pattern)) if (text[i + j] & pattern[j]) == 0) if modifier else sum(1 for j in range(len(pattern)) if text[i + j] != pattern[j])
        if mismatch <= tolerance:
            score = (len(pattern) - mismatch) / len(pattern)
            matches.add((f'{i + 1}..{i + len(pattern)}', score, mismatch))
    return matches

def map_short_to_long(short: Seq, long: Seq, sim_thres: float, cov_thres: float, is_nucl: bool, is_circ: bool) -> list[tuple] | None:
    if is_circ and short.length() * 2 > long.length(): short, long = long, short
    elif short.length() > long.length(): short, long = long, short
    max_mismatch = short.length() - round(sim_thres * short.length())
    coverage = short.length() / long.length()
    if coverage < cov_thres: return None
    results = set()
    def process_matches(matches, strand, sequence, other_matches=None):
        for region, mm, score in matches:
            start, end = map(int, region.split('..'))
            if end > long.length() and start < long.length(): end -= long.length(); region = f'{start}..{end} (circular)'
            elif start > long.length(): continue
            match = (short.id, long.id, region, mm, score, coverage, strand)
            if other_matches:
                for o_region, o_mm, o_score in other_matches:
                    o_start, o_end = map(int, o_region.split('..'))
                    if o_end > long.length() and o_start < long.length(): o_end -= long.length(); o_region = f'{o_start}..{o_end} (circular)'
                    elif o_start > long.length(): continue
                    if o_score > score and o_mm < mm: match = (short.id, long.id, o_region, o_mm, o_score, coverage, '-' if strand == '+' else '+')
            results.add(match)

    if is_nucl:
        rc_short = reverse_complement(short.sequence)
        long_sequence = long.sequence if not is_circ else long.sequence + long.sequence
        fw_matches = search(long_sequence, short.sequence, dna_to_bits, max_mismatch)
        rc_matches = search(long_sequence, rc_short, dna_to_bits, max_mismatch)
        if fw_matches: process_matches(fw_matches, '+', short.sequence, rc_matches)
        elif rc_matches: process_matches(rc_matches, '-', rc_short)
    else:
        matches = search(long.sequence, short.sequence, None, max_mismatch)
        process_matches(matches, '.', short.sequence)
    return [(match[0], match[1], match[2], f'{match[3] * 100:.2f}', match[4], f'{match[5] * 100:.2f}', match[6]) for match in results] if results else None

def map_sequences(query_sequences: list[Seq], reference_sequences: list[Seq], similarity_threshold: float, coverage_threshold: float, is_nucleotide: bool, is_circular: bool, output_file: str, num_threads: int) -> None:
    chunk_size = 100
    results = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        func = partial(map_short_to_long, sim_thres=similarity_threshold, cov_thres=coverage_threshold, is_nucl=is_nucleotide, is_circ=is_circular)
        for chunk in chunk_iterable(product(query_sequences, reference_sequences), chunk_size):
            chunk_futures = [executor.submit(func, query, reference) for query, reference in chunk]
            futures.extend(chunk_futures)
        for future in as_completed(futures):
            if future.result(): results.extend(future.result())
    headers = ['Query_ID', 'Reference_ID', 'Region', 'Mismatch', 'Similarity', 'Coverage', 'Strand']
    if results:
        results = sorted(results, key=lambda x: x[4], reverse=True)
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)
            csv_writer.writerows(results)
    else:
        logging.info(f'No matches found.')

def chunk_iterable(iterable: list, chunk_size: int):
    iterator = iter(iterable)
    for first in iterator:  # Start each chunk with the first item
        chunk = [first] + list(islice(iterator, chunk_size - 1))
        if not chunk: break
        yield chunk

def main():
    parser = argparse.ArgumentParser(description='Simple sequence mapping.')
    parser.add_argument('--query', required=True, help='Path to the query input file (FASTA/FASTQ)')
    parser.add_argument('--reference', required=True, help='Path to the reference input file (FASTA/FASTQ)')
    parser.add_argument('-s', '--similarity', type=float, default=0.8, help='Similarity threshold for sequence matching (default: 0.8)')
    parser.add_argument('-c', '--coverage', type=float, default=0.0, help='Coverage threshold for sequence matching (default: 0.0)')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file')
    parser.add_argument('-m', '--mode', choices=['nu', 'aa'], default='nu', help='Comparison mode: "nu" for DNA/RNA, "aa" for proteins (default: "nu")')
    parser.add_argument('--circular', action='store_true', help='Set if the reference sequence is circular')
    args = parser.parse_args()
    query_sequences = read_sequences(args.query)
    reference_sequences = read_sequences(args.reference)
    is_nucleotide = args.mode == 'nu'
    is_circular = args.circular
    map_sequences(query_sequences, reference_sequences, args.similarity, args.coverage, is_nucleotide, is_circular, args.output, args.threads)

if __name__ == '__main__':
    main()
