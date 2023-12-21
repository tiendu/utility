import resource
import os
import sys
import concurrent.futures
import gzip
import zipfile
import re
import argparse
from functools import partial, lru_cache
from dataclasses import dataclass
from itertools import product, groupby
import logging
from typing import List, Tuple, Any

# Setup logging to display messages with INFO level and above
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __hash__(self):
        return hash((self.id, self.sequence))

def gen_regex(primer: str) -> str:
    nuc_variants = {
        'A': ['A'],
        'T': ['T'],
        'G': ['G'],
        'C': ['C'],
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
        'N': ['A', 'T', 'G', 'C']
    }
    primer_bases = list(primer)
    regex_pattern = ''
    for base in primer_bases:
        variants = nuc_variants.get(base, [base])
        regex_pattern += '[' + ''.join(variants) + ']' if len(variants) > 1 else re.escape(base)
    return regex_pattern

"""
@lru_cache(maxsize=None)
def gen_variants(seq: str) -> List[str]:
    nuc_vars = {
        'A': ['A'],
        'T': ['T'],
        'G': ['G'],
        'C': ['C'],
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
        'N': ['A', 'T', 'G', 'C']
    }
    var_lists = [nuc_vars.get(nuc, [nuc]) for nuc in seq]
    vars_gen = [''.join(var) for var in product(*var_lists)]
    return vars_gen
"""

def rev_complement(dna_seq: str) -> str:
    comp_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp_seq = ''.join(comp_bases.get(base, base) for base in reversed(dna_seq))
    return rev_comp_seq

def check_primers(sequence: Seq, fwd_primer_regex: Any, rev_primer_regex: Any) -> Tuple[str, str, str]:
    fwd_match = fwd_primer_regex.search(sequence.sequence)
    rev_match = rev_primer_regex.search(rev_complement(sequence.sequence))
    if fwd_match and rev_match:
        return sequence.id, fwd_match.group(), rev_match.group()
    return None, None, None

def check_primers_batch(sequences: List[Seq], fwd_primer_regex: Any, rev_primer_regex: Any) -> List[Tuple[str, str, str]]:
    return [check_primers(sequence, fwd_primer_regex, rev_primer_regex) for sequence in sequences]

def read_sequences_from_file(file_path: str, file_type: str) -> List[Seq]:
    sequences = []

    # Determine the appropriate file opening mode based on the file extension.
    opener = gzip.open if file_path.endswith('.gz') else open

    if file_type == "FASTQ":
        with opener(file_path, 'rt') as fin:
            groups = groupby(enumerate(fin), key=lambda x: x[0] // 4)
            for _, group in groups:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                name = header_line[1:]  # Removing "@" character
                seq = sequence_line.upper()
                qual = quality_line.upper()
                sequences.append(Seq(name, seq, qual))
    elif file_type == "FASTA":
        with opener(file_path, 'rt') as fin:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == ">"))
            for header in faiter:
                headerstr = next(header).strip()
                name = headerstr[1:]  # Removing ">" character
                seq = "".join(s.strip().upper() for s in next(faiter))
                sequences.append(Seq(name, seq))
    return sequences
    
def write_sequences_to_file(sequences: List[Seq], file_path: str) -> None:
    with open(file_path, 'w') as f:
        for seq in sequences:
            if seq.quality == '':
                f.write(f">{seq.id}\n{seq.sequence}\n")
            else:
                f.write(f"@{seq.id}\n{seq.sequence}\n+\n{seq.quality}\n")

def process_file(input_file: str, output_file: str, file_type: str, fwd_primer: str, rev_primer: str, num_threads: int) -> List[Tuple[str, str, str]]:
    logging.info(f"Processing file: {input_file}")
    try:
        sequences = read_sequences_from_file(input_file, file_type)
    except Exception as e:
        raise ValueError(f'An error occurred while reading the file: {str(e)}')
        
    logging.info("Compiling regex patterns.")
    fwd_primer_regex = re.compile(gen_regex(fwd_primer))
    rev_primer_regex = re.compile(gen_regex(rev_primer))
    logging.info("Regex patterns compiled.")
    
    total_sequences = len(sequences)
    chunk_size = max(1, total_sequences // num_threads)
    seq_chunks = [sequences[i:i + chunk_size] for i in range(0, total_sequences, chunk_size)]
    
    check_batch_with_primers = partial(check_primers_batch, fwd_primer_regex=fwd_primer_regex, rev_primer_regex=rev_primer_regex)
    logging.info("Checking for primer sequences.")
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = executor.map(check_batch_with_primers, seq_chunks)
        matched_ids_primers = [res for batch in results for res in batch if res[0] is not None]
    matched_seqs = []
    for sequence in sequences:
        for match in matched_ids_primers:
            if sequence.id == match[0]: 
                matched_seqs.append(Seq(f"{sequence.id}_{match[1]}_{match[2]}", sequence.sequence, sequence.quality))
    write_sequences_to_file(matched_seqs, output_file)
    	    
def main():
    parser = argparse.ArgumentParser(description='Process DNA sequences with primers.')
    parser.add_argument('-i', '--input_file', required=True, type=str, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, type=str, help='Path to the output file.')
    parser.add_argument('-f', '--forward_primer', required=True, type=str, help='Forward primer sequence')
    parser.add_argument('-r', '--reverse_primer', required=True, type=str, help='Reverse primer sequence')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads (default: 4)')
    args = parser.parse_args()

    file_type = ''
    if any(ext in args.input_file for ext in ['.fastq', '.fq']):
        file_type = "FASTQ"
    elif any(ext in args.input_file for ext in ['.fasta', '.fa', '.fna']):
        file_type = "FASTA"
    else:
        raise ValueError(f"Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).")
        
    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}.")
        args.num_threads = min(max(args.num_threads, 1), max_threads)
    if args.num_threads >= len(read_sequences_from_file(args.input_file, file_type)) * 0.1:
        logging.warning(f'Number of sequences too low. Adjusting thread count to 1.')
        args.num_threads = 1

    try:
        process_file(args.input_file, args.output_file, file_type, args.forward_primer, args.reverse_primer, args.num_threads)
    except Exception as e:
        logging.error(f"An error occurred during processing: {str(e)}")

if __name__ == '__main__':
    main()
