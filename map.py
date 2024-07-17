from typing import List, Dict
from dataclasses import dataclass
import concurrent.futures
import os
import re
from functools import partial
from itertools import product, groupby
import sys
import gzip

@dataclass(frozen=False)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __hash__(self):
        return hash((self.id, self.sequence, self.quality))

FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']
DISPLAY_LENGTH = 20

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

def reverse_translate(sequence: str) -> str:
    conversion = {
        'W': 'TGG', 'Y': 'TAY', 'C': 'TGY', 'E': 'GAR',
        'K': 'AAR', 'Q': 'CAR', 'S': 'WSN', 'L': 'YTN',
        'R': 'MGN', 'G': 'GGN', 'F': 'TTY', 'D': 'GAY',
        'H': 'CAY', 'N': 'AAY', 'M': 'ATG', 'A': 'GCN',
        'P': 'CCN', 'T': 'ACN', 'V': 'GTN', 'I': 'ATH',
        '*': 'TRR', 'X': 'NNN'
    }
    converted = [conversion[aa] for aa in sequence]
    return ''.join(converted)

def create_regex(sequence: str) -> str:
    ambiguous_nucleotide_patterns = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    regex_pattern = []
    try:
        for nu in sequence:
            regex_pattern.append(ambiguous_nucleotide_patterns[nu])
    except KeyError as e:
        raise ValueError(f"Invalid nucleotide '{e.args[0]}' in sequence. Valid nucleotides: {list(ambiguous_nucleotide_patterns.keys())}") from e
    return ''.join(regex_pattern)

def index_kmers(string: str, k: int) -> Dict[str, List[int]]:
    '''
    Index a string with k-mers of length k,
    returning a dictionary where keys are k-mers and values
    are lists of indices.
    '''
    kmer_dict = {}
    for i in range(len(string) - k + 1):
        kmer = string[i:i+k]
        if kmer in kmer_dict:
            kmer_dict[kmer].append(i)
        else:
            kmer_dict[kmer] = [i]
    return kmer_dict

def truncate_string(string: str, max_length: int) -> str:
    if len(string) > max_length:
        return string[:max_length-3] + '...'
    return string

def match_pair(sequence1: Seq, sequence2: Seq, sequence1_type: str, sequence2_type: str) -> List:
    matches = []

    if sequence1_type == 'aa':
        sequence1.sequence = reverse_translate(sequence1.sequence)
    elif sequence1_type == 'nu':
        pass
    else:
        raise ValueError('Incorrect input 1 type!')

    if sequence2_type == 'aa':
        sequence2.sequence = reverse_translate(sequence2.sequence)
    elif sequence2_type == 'nu':
        pass
    else:
        raise ValueError('Incorrect input 2 type!')

    if len(sequence1.sequence) > len(sequence2.sequence):
        k = len(sequence2.sequence)
        query = sequence2
        reference = sequence1
        reference_type = sequence1_type
    else:
        k = len(sequence1.sequence)
        query = sequence1
        reference = sequence2
        reference_type = sequence2_type

    reference_kmers = index_kmers(reference.sequence, k)
    for kmer in reference_kmers:
        if re.fullmatch(create_regex(kmer), query.sequence):
            for index in reference_kmers[kmer]:
                # print(index, kmer, query.sequence)
                start = index + 1
                end = index + k
                if reference_type == 'nu':
                    pass
                elif reference_type == 'aa':
                    start = int(start / 3) + 1
                    end = int(end / 3)
                matches.append((query.id, reference.id, f'{start}..{end}'))

    return matches

def process_concurrently(sequences1: List[Seq], sequences2: List[Seq], sequences1_type: str, sequences2_type: str) -> None:
    results = []

    max_workers = os.cpu_count()

    # Use ProcessPoolExecutor to parallelize the execution
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
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

        # Determine the maximum width for each column
        max_widths = [max(len(header), max(len(row[i]) for row in truncated_results)) for i, header in enumerate(headers)]

        # Print table headers
        print('|'.join(f'{header.ljust(max_widths[i])}' for i, header in enumerate(headers)))
        print('|'.join('-' * width for width in max_widths))

        # Print table rows
        for row in truncated_results:
            print('|'.join(f'{cell.ljust(max_widths[i])}' for i, cell in enumerate(row)))

def main():
    if len(sys.argv) < 5:
        print(f'Usage: python {sys.argv[0]} <input_1> <input_2> <input_1_type> <input_2_type>')
        return

    file1, file2, file1_type, file2_type = sys.argv[1:5]
    file1 = read_sequences(file1)
    file2 = read_sequences(file2)
    process_concurrently(file1, file2, file1_type, file2_type)

if __name__ == '__main__':
    main()
