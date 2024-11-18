import argparse
import os
import gzip
import logging
import hashlib
from itertools import groupby
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']

# Setup logging to display messages with INFO level and above.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __post_init__(self):
        object.__setattr__(self, 'sequence', self.sequence.upper())

    def __hash__(self) -> int:
        return hash((self.id, self.sequence, self.quality))

    def length(self) -> int:
        return len(self.sequence)

def hash_string(string: str, hash_function=hashlib.md5) -> str:
    return hash_function(string.encode()).hexdigest()

def read_sequences_from_file(file_path: str, file_type: str) -> list[Seq]:
    unique_seqs = {}
    opener = gzip.open if file_path.endswith('.gz') else open
    count = 0
    file_handler = opener(file_path, 'rt')
    if file_type == 'FASTQ':
        groups = groupby(enumerate(file_handler), key=lambda x: x[0] // 4)
        for _, group in groups:
            try:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                seqid = header_line[1:]
                seq = sequence_line
                qual = quality_line
                count += 1
                seq_hash = hash_string(seq)
                unique_seqs.setdefault(seq_hash, Seq(seqid, seq, qual))
            except ValueError:
                logging.warning('Malformed FASTQ entry encountered. Skipping')
    elif file_type == 'FASTA':
        faiter = (x[1] for x in groupby(file_handler, lambda line: line[0] == '>'))
        for header in faiter:
            try:
                header_str = next(header).strip()
                seqid = header_str[1:]
                seq = ''.join(s.strip() for s in next(faiter))
                count += 1
                seq_hash = hash_string(seq)
                unique_seqs.setdefault(seq_hash, Seq(seqid, seq))
            except StopIteration:
                logging.warning('Incomplete FASTA entry encountered. Skipping')
    file_handler.close()
    logging.info(f'Start: {count}')
    return sorted(unique_seqs.values(), key=lambda x: x.length())

def write_sequences_to_file(sequences: list[Seq], file_path: str) -> None:
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'wt') as f:
        for seq in sequences:
            f.write(f'{">" if seq.quality == "" else "@"}{seq.id}\n{seq.sequence}\n')
            if seq.quality:
                f.write(f'+\n{seq.quality}\n')

def deduplicate_sequences(seqs: list[Seq], num_threads: int) -> list[Seq]:
    length_change_indices = [i for i, seq in enumerate(seqs) if seq.length() > (seqs[i - 1].length() if i > 0 else -1)]
    total_sequences = len(seqs)
    progress_interval = max(1, total_sequences // 20)
    def process_sequence(short, i):
        next_index = next((index for index in length_change_indices if index > i), total_sequences)
        seqs[i] = None if any(short.sequence in longer.sequence for longer in seqs[next_index:]) else seqs[i]

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = {executor.submit(process_sequence, seq, i): i for i, seq in enumerate(seqs)}
        for progress_counter, _ in enumerate(as_completed(futures), start=1):
            if progress_counter % progress_interval == 0:
                logging.info(f'Progress: {progress_counter / total_sequences:.1%} ({progress_counter}/{total_sequences})')
    seqs = [seq for seq in seqs if seq is not None]
    logging.info(f'End: {len(seqs)}')
    return seqs

def main():
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')
    args = parser.parse_args()
    if not os.path.exists(args.input_file):
        raise ValueError(f'The input file "{args.input_file}" does not exist.')
    try:
        with open(args.output_file, 'w'):
            pass
    except Exception as e:
        raise ValueError(f'Error creating output file "{args.output_file}": {e}')
    file_type = (
        'FASTQ' if any(ext in args.input_file for ext in FASTQ_EXTENSIONS)
        else 'FASTA' if any(ext in args.input_file for ext in FASTA_EXTENSIONS)
        else ValueError(f'Unrecognized file extension for {args.input_file}.')
    )
    args.num_threads = (
        max(1, min(args.num_threads, os.cpu_count()))
        if args.num_threads > 0
        else os.cpu_count()
    )
    sequences = read_sequences_from_file(args.input_file, file_type)
    if not sequences:
        raise ValueError('No sequences detected!')
    args.num_threads = max(1, min(args.num_threads, len(sequences) // 10))
    deduplicated_sequences = deduplicate_sequences(sequences, args.num_threads)
    write_sequences_to_file(deduplicated_sequences, args.output_file)

if __name__ == '__main__':
    main()
