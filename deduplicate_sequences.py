import argparse
import os
import gzip
import logging
import hashlib
from itertools import groupby
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor

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
    def __hash__(self) -> int:
        return hash((self.id, self.sequence, self.quality))

    def length(self) -> int:
        return len(self.sequence)

def hash_string(string: str, hash_function=hashlib.md5) -> str:
    return hash_function(string.encode()).hexdigest()

def add_sequence_to_collections(uniq_seqs: dict, sorted_seqs: list,
                                identifier: str, sequence: str,
                                quality: str = '') -> None:
    def insert_in_sorted_order(sorted_seqs: list[Seq], new_seq: Seq) -> None:
        new_seq_len = new_seq.length()
        for i, existing_seq in enumerate(sorted_seqs):
            if existing_seq.length() > new_seq_len:
                sorted_seqs.insert(i, new_seq)
                return
        sorted_seqs.append(new_seq)

    seq_hash = hash_string(sequence)
    if seq_hash not in uniq_seqs:
        uniq_seqs[seq_hash] = Seq(identifier, sequence, quality)
        insert_in_sorted_order(sorted_seqs, uniq_seqs[seq_hash])

def read_sequences_from_file(file_path: str, file_type: str) -> list[Seq]:
    uniq_seqs = {}
    sorted_seqs = []
    opener = gzip.open if file_path.endswith('.gz') else open
    if file_type == 'FASTQ':
        with opener(file_path, 'rt') as fin:
            groups = groupby(enumerate(fin), key=lambda x: x[0] // 4)
            for _, group in groups:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                seqid = header_line[1:]
                seq = sequence_line.upper()
                qual = quality_line
                add_sequence_to_collections(uniq_seqs, sorted_seqs, seqid, seq, qual)
    elif file_type == 'FASTA':
        with opener(file_path, 'rt') as fin:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == '>'))
            for header in faiter:
                header_str = next(header).strip()
                seqid = header_str[1:]
                seq = ''.join(s.strip().upper() for s in next(faiter))
                add_sequence_to_collections(uniq_seqs, sorted_seqs, seqid, seq)
    return sorted_seqs

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

def hash_string(string: str, hash_function=hashlib.md5) -> str:
    return hash_function(string.encode()).hexdigest()

def dna_to_bits(dna: str) -> list[int]:
    nu_to_bits = {
        'A': 0b0001, 'T': 0b0010, 'G': 0b0100, 'C': 0b1000,
        'W': 0b0011, 'R': 0b0101, 'Y': 0b1010, 'S': 0b1100,
        'K': 0b0110, 'M': 0b1001, 'B': 0b1110, 'D': 0b0111,
        'H': 0b1011, 'V': 0b1101, 'N': 0b1111
    }
    return [nu_to_bits[nu.upper()] for nu in dna]

def compute_lps(pattern: str) -> list[int]:
    lps = [0] * len(pattern)
    length = 0  # Length of the previous longest prefix suffix
    i = 1
    while i < len(pattern):
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1
    return lps

def kmp_search(text: str, pattern: str, lps: list[int], modifier=None) -> bool:
    def match(a: str | int, b: str | int) -> bool:
        if isinstance(a, int) and isinstance(b, int):
            return (a & b) != 0
        return a == b

    text = modifier(text) if modifier else text
    pattern = tuple(modifier(pattern)) if modifier else pattern
    i = 0  # Index for text
    j = 0  # Index for pattern
    while i < len(text):
        if match(pattern[j], text[i]):
            i += 1
            j += 1
        if j == len(pattern):
            return True  # Pattern found
        elif i < len(text) and not match(pattern[j], text[i]):
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return False  # Pattern not found

def check_sequence(short: Seq, longers: list[Seq], is_nucl: bool) -> Seq | None:
    lps = compute_lps(short.sequence)  # Compute LPS array
    modifier = dna_to_bits if is_nucl else None
    for longer in longers:
        if kmp_search(longer.sequence, short.sequence, lps, modifier):
            return None
    return short

def deduplicate_sequences(sequences: list[Seq], num_threads: int, is_nucl: bool) -> list[Seq]:
    short_to_long = sorted(sequences,
                           key=lambda seq: seq.length(),
                           reverse=False)
    deduped_set = set()
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for i, short in enumerate(short_to_long):
            longers = [seq for seq in short_to_long[i + 1:]
                       if seq.length() > short.length()]
            future = executor.submit(check_sequence,
                                     short,
                                     longers,
                                     is_nucl)
            futures.append(future)
        for future in futures:
            result = future.result()
            if result:
                deduped_set.add(result)
                logging.info(f'No. of deduped sequences: {len(deduped_set)}')
    return deduped_set

def main() -> None:
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')
    parser.add_argument('-m', '--mode', choices=['nu', 'aa'], default='nu',
                        help="Comparison mode: 'nu' for DNA/RNA, 'aa' for proteins (default: 'nu')")
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
    is_nucleotide = True if args.mode == 'nu' else False
    deduplicated_sequences = deduplicate_sequences(sequences, args.num_threads, is_nucleotide)
    write_sequences_to_file(deduplicated_sequences, args.output_file)

if __name__ == '__main__':
    main()
