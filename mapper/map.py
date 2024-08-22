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

def get_minimizer(sequence: str, k: int, w: int) -> List[Tuple[str, int]]:
    minimizers = []
    n = len(sequence)
    
    for i in range(n - w + 1):
        window = sequence[i:i + w]
        k_mers = [(window[j:j + k], i + j) for j in range(w - k + 1)]
        min_kmer = min(k_mers)
        minimizers.append(min_kmer)

    unique_minimizers = []
    for m in minimizers:
        if not unique_minimizers or m[0] != unique_minimizers[-1][0]:
            unique_minimizers.append(m)
    
    return unique_minimizers

def match_minimizers(minimizers1: List[Tuple[str, int]], minimizers2: List[Tuple[str, int]], threshold: float) -> float:
    set1 = set(minimizer[0] for minimizer in minimizers1)
    set2 = set(minimizer[0] for minimizer in minimizers2)
    
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    
    if not union:
        return 0.0
    
    similarity = len(intersection) / len(union)
    
    if similarity > threshold:
        return similarity

def map(query: Dict[str, Dict[str, str]], database: List[Seq], k: int, w: int, threshold: float) -> Dict[str, List[Tuple[str, float]]]:
    matches = {}

    for db_seq in database:
        for index, seq_info in query.items():
            k = len(index)
            for i in range(len(sequence) - k + 1):
                db_kmer = sequence[i:i + k]
                db_kmer_minimizers = get_minimizer(kmer, len(index) // 2 | 1, len(index))
                index_minimizers = get_minimizer(index, len(index) // 2 | 1, len(index))
                index_similarity = match_minimizers(index_minimizers, db_kmer_minimizers, threshold)
                if index_similarity == 0:
                    continue
            db_minimizers = get_minimizer(db_seq.sequence, len(index) // 2 | 1, len(index))
            query_minimizers = get_minimizer(seq.sequence, len(index) // 2 | 1, len(index))
            similarity = match_minimizers(query_minimizers, db_minimizers, threshold)
                    
    return matches

def main():
    parser = argparse.ArgumentParser(description='Process sequences from FASTA/FASTQ files.')
    parser.add_argument('-q', '--query', type=str, required=True, help='Path to the indexed query file (JSON).')
    parser.add_argument('-d', '--database', type=str, required=True, help='Path to the database file (FASTA/FASTQ).')
    parser.add_argument('-m', '--mode', choices=['aa', 'nu'], required=True, help="Mode: 'aa' for amino acid sequences, 'nu' for nucleotide sequences.")
    parser.add_argument('-o', '--output', type=str, help='Output file path. Defaults to <input_file>.json')
    args = parser.parse_args()

    query = Path(args.query)
    with query.open('r') as f:
        query = json.load(f)
    
    outfile = Path(args.output) if args.output else infile.with_suffix('.json')

    seqs = read_sequences(infile)

    for seq in seqs:
        if args.mode == 'aa':
            seq.sequence = reverse_translate(seq.sequence)

    delineated_sequences = []
    with ProcessPoolExecutor() as executor:
        future_to_seq = {executor.submit(process_sequence, seq, args.mode): seq for seq in seqs}
        for future in as_completed(future_to_seq):
            result = future.result()
            delineated_sequences.extend(result)

    k = min(len(seq.sequence) // 2 | 1 for seq in delineated_sequences)
    kmer_index = index(delineated_sequences, k)

    with outfile.open('w') as outfile:
        json.dump({kmer: {'id': seq.id, 'sequence': seq.sequence, 'quality': seq.quality} 
                   for kmer, seq in kmer_index.items()}, outfile, indent=4)

if __name__ == '__main__':
    main()
