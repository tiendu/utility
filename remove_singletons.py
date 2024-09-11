import argparse
from dataclasses import dataclass
from collections import Counter
from pathlib import Path
from typing import Generator, List
import gzip
import sys

# Constants
FASTQ_EXTENSIONS = {'.fastq', '.fq'}
FASTA_EXTENSIONS = {'.fasta', '.fa', '.fna', '.faa'}
LINE_LENGTH = 80  # Standard line length for FASTA formatting

@dataclass(frozen=False)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    
    def __post_init__(self):
        object.__setattr__(self, 'sequence', self.sequence.upper())

def read_sequences(file_path: Path) -> List[Seq]:
    sequences = []
    file_type = None
    if file_path.suffix in FASTQ_EXTENSIONS:
        file_type = 'FASTQ'
    elif file_path.suffix in FASTA_EXTENSIONS:
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {file_path}. Expected FASTA {FASTA_EXTENSIONS} or FASTQ {FASTQ_EXTENSIONS}')
    
    opener = gzip.open if file_path.suffix == '.gz' else open
    with opener(file_path, 'rt') as fin:
        if file_type == 'FASTQ':
            while True:
                header_line = fin.readline().strip()
                if not header_line:
                    break
                sequence_line = fin.readline().strip()
                _ = fin.readline().strip()  # Plus line
                quality_line = fin.readline().strip()
                sequences.append(Seq(header_line[1:], sequence_line, quality_line))
        elif file_type == 'FASTA':
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

def write_sequences_to_file(sequences: List[Seq], file_path: Path) -> None:
    opener = gzip.open if file_path.suffix == '.gz' else open
    with opener(file_path, 'wt') as f:
        for seq in sequences:
            if seq.quality == '':
                f.write(f'>{seq.id}\n')
                for i in range(0, len(seq.sequence), LINE_LENGTH):
                    f.write(seq.sequence[i:i+LINE_LENGTH] + '\n')
            else:
                f.write(f'@{seq.id}\n{seq.sequence}\n+\n{seq.quality}\n')

def kmerize(string: str, k: int) -> List[str]:
    return [string[i:i + k] for i in range(len(string) - k + 1)]

def remove_singletons(sequences: List[Seq]) -> List[Seq]:
    min_length = min(len(sequence.sequence) for sequence in sequences)
    
    kmer_count = Counter()
    sequence_kmers = {}
    
    for seq in sequences:
        kmers = kmerize(seq.sequence, min_length)
        sequence_kmers[seq.id] = kmers
        kmer_count.update(kmers)
    
    filtered_sequences = []
    for seq in sequences:
        if any(kmer_count[kmer] > 1 for kmer in sequence_kmers[seq.id]):
            filtered_sequences.append(seq)
    
    return filtered_sequences

def main():
    parser = argparse.ArgumentParser(description='Remove singletons using k-mers.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file')
    args = parser.parse_args()
    
    input_path = Path(args.input)
    output_path = Path(args.output)
    
    sequences = read_sequences(input_path)
    filtered_sequences = remove_singletons(sequences)
    
    write_sequences_to_file(filtered_sequences, output_path)

if __name__ == '__main__':
    main()
