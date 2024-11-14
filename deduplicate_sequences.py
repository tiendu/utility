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
LINE_LENGTH = 80

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
    if file_type == 'FASTQ':
        with opener(file_path, 'rt') as fin:
            groups = groupby(enumerate(fin), key=lambda x: x[0] // 4)
            for _, group in groups:
                try:
                    header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                    seqid = header_line[1:]
                    seq = sequence_line
                    qual = quality_line
                    seq_hash = hash_string(seq)
                    if seq_hash not in unique_seqs:
                        count += 1
                        unique_seqs[seq_hash] = Seq(seqid, seq, qual)
                except ValueError:
                    logging.warning('Malformed FASTQ entry encountered. Skipping')
    elif file_type == 'FASTA':
        with opener(file_path, 'rt') as fin:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == '>'))
            for header in faiter:
                try:
                    header_str = next(header).strip()
                    seqid = header_str[1:]
                    seq = ''.join(s.strip() for s in next(faiter))
                    seq_hash = hash_string(seq)
                    if seq_hash not in unique_seqs:
                        count += 1
                        unique_seqs[seq_hash] = Seq(seqid, seq)
                except StopIteration:
                    logging.warning('Incomplete FASTA entry encountered. Skipping')
    logging.info(f'Read sequences: {count}')
    return sorted(unique_seqs.values(), key=lambda x: x.length())

def write_sequences_to_file(sequences: list[Seq], file_path: str) -> None:
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'wt') as f:
        for seq in sequences:
            if seq.quality == '':
                f.write(f'>{seq.id}\n')
                for i in range(0, len(seq.sequence), LINE_LENGTH):
                    f.write(seq.sequence[i:i+LINE_LENGTH] + '\n')
            else:
                f.write(f'@{seq.id}\n{seq.sequence}\n+\n{seq.quality}\n')

def check_sequence(short: Seq, longers: list[Seq]) -> Seq | None:
    for longer in longers:
        if short.sequence in longer.sequence:
            return None
    return short
def deduplicate_sequences(seqs: list[Seq], num_threads: int) -> list[Seq]:
    deduped = set()
    length_change_indices = []
    last_length = -1
    for i, seq in enumerate(seqs):  # Calculate indices where sequence length changes
        if seq.length() > last_length:
            length_change_indices.append(i)
            last_length = seq.length()
    total_sequences = len(seqs)
    progress_interval = max(1, total_sequences // 20)  # Log every 5%
    progress_counter = 0  # Tracks completed futures for logging progress
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for i, short in enumerate(seqs):
            next_index = next((index for index in length_change_indices if index > i), total_sequences)  # Identify the range of longer sequences
            longers = seqs[next_index:]
            future = executor.submit(check_sequence, short, longers)
            futures.append(future)
        for future in as_completed(futures):
            result = future.result()
            if result:
                deduped.add(result)
            progress_counter += 1
            if progress_counter % progress_interval == 0:
                progress_percent = (progress_counter / total_sequences) * 100
                logging.info(f'Progress: {progress_percent:.1f}% ({progress_counter}/{total_sequences})')
    logging.info(f'Total deduplicated sequences: {len(deduped)}')
    return deduped

def main() -> None:
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
    file_type = ''
    if any(ext in args.input_file for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(ext in args.input_file for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna, .faa) or FASTQ (.fastq, .fq).')
    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f'Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}')
        args.num_threads = min(max(args.num_threads, 1), max_threads)
    sequences = read_sequences_from_file(args.input_file, file_type)
    if not sequences:
        raise ValueError('No sequences detected!')
    if args.num_threads >= len(sequences) * 0.1:
        logging.warning(f'Number of sequences too low, thread count adjusted to 1')
        args.num_threads = 1
    deduplicated_sequences = deduplicate_sequences(sequences, args.num_threads)
    write_sequences_to_file(deduplicated_sequences, args.output_file)

if __name__ == '__main__':
    main()
