import sys
import re
import os
import gzip
from itertools import groupby, product
from dataclasses import dataclass
from typing import List
from concurrent.futures import ProcessPoolExecutor
import hashlib
from concurrent.futures import as_completed
import numpy as np

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __hash__(self):
        return hash((self.id, self.sequence, self.quality))

FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']

def read_sequences(file_path: str) -> List[Seq]:
    sequences = []

    file_type = ''
    if any(ext in file_path for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(ext in file_path for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {file_path}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).')
    
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

def generate_kmers(string: str, k: int) -> np.ndarray:
    kmers = [string[i:i+k] for i in range(len(string) - k + 1)]
    return np.array(kmers)

def hash_sequence(sequence: str, hash_function=hashlib.sha3_256) -> str:
    return hash_function(sequence.encode()).hexdigest()

def deduplicate(sequences: List[Seq]) -> List[Seq]:
    sequences.sort(key=lambda s: len(s.sequence), reverse=True)
    unique_kmer_hashes = set()
    unique_sequences = dict()

    if sequences:
        min_length = len(sequences[-1].sequence)

    for current_sequence in sequences:
        kmers = generate_kmers(current_sequence.sequence, min_length)
        kmer_hashes = {hash_sequence(kmer) for kmer in kmers}
        if all(hash_value in unique_kmer_hashes for hash_value in kmer_hashes):
            continue

        unique_sequences[hash_sequence(current_sequence.sequence)] = current_sequence
        unique_kmer_hashes.update(kmer_hashes)

    return list(unique_sequences.values())

def reverse_translate(sequences: List[Seq]) -> List[Seq]:
    conversion = {
        'W': 'TGG', 'Y': 'TAY', 'C': 'TGY', 'E': 'GAR',
        'K': 'AAR', 'Q': 'CAR', 'S': 'WSN', 'L': 'YTN',
        'R': 'MGN', 'G': 'GGN', 'F': 'TTY', 'D': 'GAY',
        'H': 'CAY', 'N': 'AAY', 'M': 'ATG', 'A': 'GCN',
        'P': 'CCN', 'T': 'ACN', 'V': 'GTN', 'I': 'ATH',
        '*': 'TRR', 'X': 'NNN'
    }
    reverse_translated = []
    for sequence in sequences:
        converted = [conversion[aa] for aa in sequence.sequence]
        reverse_translated.append(Seq(id=sequence.id, sequence=''.join(converted)))
        
    return reverse_translated

def create_regex(sequence: str) -> re.Pattern:
    ambiguous_nucleotide_patterns = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    regex_pattern = ''.join(ambiguous_nucleotide_patterns[nu] for nu in sequence)
    return re.compile(regex_pattern)

def match_pair(sequence1: Seq, sequence2: Seq) -> List:
    matches = []
    
    if len(sequence1.sequence) > len(sequence2.sequence):
        min_length = len(sequence2.sequence)
        max_length = len(sequence1.sequence)
        query = sequence2
        reference = sequence1
    else:
        min_length = len(sequence1.sequence)
        max_length = len(sequence2.sequence)
        query = sequence1
        reference = sequence2

    query_pattern = create_regex(query.sequence)
    for i in range(max_length - min_length + 1):
        sub_reference = reference.sequence[i:i + min_length]
        if len(sub_reference) == min_length:
            sub_reference_pattern = create_regex(sub_reference)
            if (query_pattern.fullmatch(sub_reference) or
                sub_reference_pattern.fullmatch(query.sequence) or
                re.fullmatch(sub_reference, query.sequence) or
                re.fullmatch(query.sequence, sub_reference)):
                matches.append((query.id, reference.id, f'{i+1}..{i+min_length}'))

    return matches

def process_concurrently(sequences1: List[Seq], sequences2: List[Seq]):
    results = []

    max_workers = os.cpu_count()

    # Use ProcessPoolExecutor to parallelize the execution
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(match_pair, sequence1, sequence2) for sequence1, sequence2 in product(sequences1, sequences2)]
        for future in as_completed(futures):
            result = future.result()
            if result:
                results.extend(result)

    def truncate_string(string: str, max_length: int) -> str:
        if len(string) > max_length:
            return string[:max_length - 3] + '...'
        return string
    
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

        # Determine the maximum width for each column
        max_widths = [max(len(header), max(len(row[i]) for row in truncated_results)) for i, header in enumerate(headers)]

        # Print table headers
        print('|'.join(f'{header.ljust(max_widths[i])}' for i, header in enumerate(headers)))
        print('|'.join('-' * width for width in max_widths))

        # Print table rows
        for row in truncated_results:
            print('|'.join(f'{cell.ljust(max_widths[i])}' for i, cell in enumerate(row)))

def map(file1: str, file2: str, file1_type: str, file2_type: str):
    sequences1 = read_sequences(file1)
    sequences2 = read_sequences(file2)

    if file1_type == 'aa':
        sequences1 = reverse_translate(deduplicate(sequences1))
    elif file1_type == 'nu':
        sequences1 = deduplicate(sequences1)
    
    if file2_type == 'aa':
        sequences2 = reverse_translate(deduplicate(sequences2))
    elif file2_type == 'nu':
        sequences2 = deduplicate(sequences2)

    process_concurrently(sequences1, sequences2)

def main():
    if len(sys.argv) < 5:
        print(f'Usage: python {sys.argv[0]} input1 input2 file1_type file2_type')
        return
    
    file1, file2, file1_type, file2_type = sys.argv[1:5]
    map(file1, file2, file1_type, file2_type)

if __name__ == "__main__":
    main()
