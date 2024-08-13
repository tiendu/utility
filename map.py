from typing import List, Dict
from dataclasses import dataclass
import concurrent.futures
import os
import re
from functools import partial
from itertools import product, groupby
import sys
import csv
import gzip

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    
    def __post_init__(self):
        object.__setattr__(self, 'sequence', self.sequence.upper())

def read_sequences(file_path: str) -> List[Seq]:
    sequences = []
    file_type = ''

    if any(ext in file_path for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(ext in file_path for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {file_path}. Expected FASTA {FASTA_EXTENSIONS} or FASTQ {FASTQ_EXTENSIONS}')

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

def reverse_translate(sequence: str) -> str:
    aa_to_nu = {
        'W': 'TGG', 'Y': 'TAY', 'C': 'TGY', 'E': 'GAR',
        'K': 'AAR', 'Q': 'CAR', 'S': 'WSN', 'L': 'YTN',
        'R': 'MGN', 'G': 'GGN', 'F': 'TTY', 'D': 'GAY',
        'H': 'CAY', 'N': 'AAY', 'M': 'ATG', 'A': 'GCN',
        'P': 'CCN', 'T': 'ACN', 'V': 'GTN', 'I': 'ATH',
        '*': 'TRR', 'X': 'NNN'
    }
    trans_table = str.maketrans(aa_to_nu)
    try:
        rev_trans = sequence.translate(trans_table)
    except KeyError as e:
        raise ValueError(f"Invalid amino acid '{e.args[0]}' in sequence. Valid amino acids: {list(aa_to_nu.keys())}") from e
    return rev_trans

def create_regex(sequence: str) -> str:
    ambig_nucl_patt = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    trans_table = str.maketrans(ambig_nucl_patt)
    try:
        regex_pattern = sequence.translate(trans_table)
    except KeyError as e:
        raise ValueError(f"Invalid nucleotide '{e.args[0]}' in sequence. Valid nucleotides: {list(ambig_nucl_patt.keys())}") from e
    return regex_pattern

def index_kmers(string: str, k: int) -> Dict[str, List[int]]:
    kmer_dict = {}
    for i in range(len(string) - k + 1):
        kmer = string[i:i+k]
        kmer_dict.setdefault(kmer, []).append(i)
    return kmer_dict

def truncate_string(string: str, max_length: int) -> str:
    return string[:max_length-3] + '...' if len(string) > max_length else string

def match_pair(sequence1: Seq, sequence2: Seq, sequence1_type: str, sequence2_type: str) -> List:
    matches = []

    if sequence1_type == 'aa':
        sequence1.sequence = reverse_translate(sequence1.sequence)
    elif sequence1_type != 'nu':
        raise ValueError('Incorrect input 1 type!')

    if sequence2_type == 'aa':
        sequence2.sequence = reverse_translate(sequence2.sequence)
    elif sequence2_type != 'nu':
        raise ValueError('Incorrect input 2 type!')

    k = min(len(sequence1.sequence), len(sequence2.sequence))
    query, reference, reference_type = (sequence2, sequence1, sequence2_type) if len(sequence1.sequence) > len(sequence2.sequence) else (sequence1, sequence2, sequence1_type)

    reference_kmers = index_kmers(reference.sequence, k)
    for kmer in reference_kmers:
        if re.fullmatch(create_regex(kmer), query.sequence):
            for index in reference_kmers[kmer]:
                start, end = index + 1, index + k
                if reference_type == 'aa':
                    start, end = (start // 3) + 1, end // 3
                matches.append((query.id, reference.id, f'{start}..{end}'))

    return matches

def process_concurrently(sequences1: List[Seq], sequences2: List[Seq], sequences1_type: str, sequences2_type: str, output_file: str) -> None:
    results = []

    max_workers = os.cpu_count() or 1  # Set the number of workers to the number of CPU cores or 1 if not available
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        func = partial(match_pair, sequence1_type=sequences1_type, sequence2_type=sequences2_type)
        futures = [executor.submit(func, sequence1, sequence2) for sequence1, sequence2 in product(sequences1, sequences2)]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.extend(result)

    if results:
        headers = ['Query ID', 'Reference ID', 'Location']
        truncated_results = [
            (
                truncate_string(query_id, 20),
                truncate_string(reference_id, 50),
                location
            )
            for query_id, reference_id, location in results
        ]

        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)
            csv_writer.writerows(truncated_results)

def main():
    if len(sys.argv) < 6:
        print(f'Usage: python {sys.argv[0]} <input_1> <input_2> <input_1_type> <input_2_type> <output_file>')
        return

    file1, file2, file1_type, file2_type, output_file = sys.argv[1:6]
    file1 = read_sequences(file1)
    file2 = read_sequences(file2)
    process_concurrently(file1, file2, file1_type, file2_type, output_file)

if __name__ == '__main__':
    main()
