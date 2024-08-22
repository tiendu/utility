from dataclasses import dataclass
import gzip
from pathlib import Path
from typing import List, Dict
from itertools import groupby
from collections import Counter
import argparse
import json

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', 'faa']

@dataclass(frozen=False)
class Seq:
    id: str
    sequence: str
    quality: str = ''

    def __post_init__(self):
        self.sequence = self.sequence.upper()

def read_sequences(file_path: Path) -> List[Seq]:
    sequences = []
    file_type = ''

    if any(file_path.suffix == ext for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(file_path.suffix == ext for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {file_path}. Expected FASTA {FASTA_EXTENSIONS} or FASTQ {FASTQ_EXTENSIONS}')

    opener = gzip.open if file_path.suffix == '.gz' else open

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
    rev_trans = []
    for aa in sequence:
        if aa in aa_to_nu:
            rev_trans.append(aa_to_nu[aa])
        else:
            raise ValueError(f"Unexpected amino acid '{aa}' encountered in sequence.")
    return ''.join(rev_trans)

def delineate(dna: str) -> List[str]:
    conversion = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'T', 'G', 'C'],
    }

    all_variants = ['']
    for nucleotide in dna:
        all_variants = [prefix + base for prefix in all_variants for base in conversion.get(nucleotide, [nucleotide])]
    
    return all_variants

def find_unique_kmer(string: str, k: int, kmer_counts: Counter) -> str:
    for i in range(len(string) - k + 1):
        kmer = string[i:i + k]
        if kmer_counts[kmer] == 1:
            return kmer
    return None

def index(sequences: List[Seq], k: int) -> Dict[str, Seq]:
    kmer_counts = Counter()

    for seq in sequences:
        sequence = seq.sequence
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            kmer_counts[kmer] += 1

    kmer_index = {}
    for seq in sequences:
        current_k = k
        unique_kmer = find_unique_kmer(seq.sequence, current_k, kmer_counts)
        while not unique_kmer and current_k < len(seq.sequence):
            current_k += 2
            unique_kmer = find_unique_kmer(seq.sequence, current_k, kmer_counts)
        
        if unique_kmer:
            kmer_index[unique_kmer] = seq
        else:
            raise ValueError(f"No unique k-mer found for sequence {seq.id}. Consider using a bigger k.")

    return kmer_index

def main():
    parser = argparse.ArgumentParser(description='Process sequences from FASTA/FASTQ files.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input sequence file (FASTA/FASTQ).')
    parser.add_argument('-m', '--mode', choices=['aa', 'nu'], required=True, help="Mode: 'aa' for amino acid sequences, 'nu' for nucleotide sequences.")
    parser.add_argument('-o', '--output', type=str, help='Output file path. Defaults to <input_file>.json')
    args = parser.parse_args()

    infile = Path(args.input)
    outfile = Path(args.output) if args.output else infile.with_suffix('.json')

    seqs = read_sequences(infile)

    for seq in seqs:
        if args.mode == 'aa':
            seq.sequence = reverse_translate(seq.sequence)

    delineated_sequences = []
    for seq in seqs:
        variants = delineate(seq.sequence)
        for i, var in enumerate(variants, 1):
            new_id = f'{seq.id}_{i}'
            delineated_sequences.append(Seq(new_id, var, seq.quality))

    k = min(len(seq.sequence) // 2 | 1 for seq in delineated_sequences)
    kmer_index = index(delineated_sequences, k)

    with outfile.open('w') as outfile:
        json.dump({kmer: {'id': seq.id, 'sequence': seq.sequence, 'quality': seq.quality} 
                   for kmer, seq in kmer_index.items()}, outfile, indent=4)

if __name__ == '__main__':
    main()
