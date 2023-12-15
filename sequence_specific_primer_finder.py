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
from itertools import product
import logging
from typing import List, Tuple, Any

logging.basicConfig(level=logging.INFO)

@dataclass
class Seq:
    id: str
    sequence: str
    quality: str = ''

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
@lru_cache(maxsize=None)
def rev_complement(dna_seq: str) -> str:
    comp_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp_seq = ''.join(comp_bases.get(base, base) for base in reversed(dna_seq))
    return rev_comp_seq

def check_primers(seq_obj: Seq, fwd_primer_regex: Any, rev_primer_regex: Any) -> Tuple[str, str, str]:
    fwd_match = fwd_primer_regex.search(seq_obj.sequence)
    rev_match = rev_primer_regex.search(rev_complement(seq_obj.sequence))
    if fwd_match and rev_match:
        return seq_obj.id, fwd_match.group(), rev_match.group()
    return None, None, None

def check_batch(seq_batch: List[Seq], fwd_primer_regex: Any, rev_primer_regex: Any) -> List[Tuple[str, str, str]]:
    return [check_primers(seq_obj, fwd_primer_regex, rev_primer_regex) for seq_obj in seq_batch]

# Function to read file
def read_file(file: str) -> List[str]:
    if file.endswith('.gz'):
        with gzip.open(file, 'rt') as gz_file:
            return gz_file.readlines()
    elif file.endswith('.qza'):
        with zipfile.ZipFile(file) as zip_file:
            sequences_entry_path = next((entry for entry in zip_file.namelist() if re.search('.*/data/dna-sequences\\.fasta', entry)), None)
            if sequences_entry_path:
                with zip_file.open(sequences_entry_path, 'r') as sequences_file:
                    return sequences_file.readlines()
            else:
                logging.error("Sequences file not found in the Qiime2 artifact.")
                raise FileNotFoundError("Sequences file not found in the Qiime2 artifact.")
    else:  # Plain text
        with open(file, 'r') as txt_file:
            return txt_file.readlines()

def process_file(input_path: str, fwd_primer: str, rev_primer: str, num_threads: int) -> List[Tuple[str, str, str]]:
    logging.info(f"Processing file: {input_path}")
    seqs_with_id = []
    if not all((isinstance(fwd_primer, str), isinstance(rev_primer, str))):
        logging.error("Wrong primers.")
        raise ValueError("All primers must be provided.")

    if not isinstance(num_threads, int) or num_threads < 1:
        logging.error("The number of threads must be a positive integer.")
        raise ValueError("The number of threads must be a positive integer.")
    # Read the file
    try:
        file_lines = read_file(input_path)
    except Exception as e:
        logging.error(f"An error occurred while reading the file: {str(e)}")
        return []
    logging.info("Compiling regex patterns.")
    fwd_primer_regex = re.compile(gen_regex(fwd_primer))
    rev_primer_regex = re.compile(gen_regex(rev_primer))
    logging.info("Regex patterns compiled.")
    is_fastq = any(input_path.endswith(extension) for extension in ['.fastq', '.fq', '.fastq.gz', '.fq.gz'])
    cur_id = None
    cur_seq = ""
    cur_quality = ""
    for line in file_lines:
        line = line.strip()
        if line.startswith(">") or line.startswith("@"):
            if cur_id and cur_seq:
                seqs_with_id.append(Seq(cur_id, cur_seq, cur_quality))
            cur_id = line[1:]
            cur_seq = ""
            cur_quality = ""
        elif is_fastq and line.startswith("+"):
            cur_quality = next(file_lines).strip()
        else:
            cur_seq += line
    if cur_id and cur_seq:
        seqs_with_id.append(Seq(cur_id, cur_seq, cur_quality))
    logging.info("Sequences extracted.")
    batch_size = -(-len(seqs_with_id) // num_threads)
    check_batch_with_primers = partial(check_batch, fwd_primer_regex=fwd_primer_regex, rev_primer_regex=rev_primer_regex)
    seq_batches = [seqs_with_id[i:i+batch_size] for i in range(0, len(seqs_with_id), batch_size)]
    logging.info("Checking for primer sequences.")
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = executor.map(check_batch_with_primers, seq_batches)
        match_ids_primers = [res for batch in results for res in batch if res[0] is not None]
    logging.info("Primer sequence check completed.")
    return match_ids_primers

def memory_limit_exceeded(limit_percent=80):
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3
    total_memory = os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE') / 1024**3
    memory_percent_used = (mem / total_memory) * 100
    if memory_percent_used > limit_percent:
        logging.error("Memory limit exceeded. Stopping execution.")
        sys.exit()

def main():
    parser = argparse.ArgumentParser(description='Process DNA sequences with primers.')
    parser.add_argument('-i', '--input_path', type=str, help='Path to the input file')
    parser.add_argument('-f', '--forward_primer', type=str, help='Forward primer sequence')
    parser.add_argument('-r', '--reverse_primer', type=str, help='Reverse primer sequence')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads (default: 4)')
    args = parser.parse_args()
    memory_limit_exceeded(80)
    try:
        matching_seqs_ids_primers = process_file(args.input_path, args.forward_primer, args.reverse_primer, args.num_threads)  # Process the file with primers
        for id, forward_primer, reverse_primer in matching_seqs_ids_primers:
            print(id, forward_primer, reverse_primer)  # Print the matching sequence ids
    except Exception as e:
        logging.error(f"An error occurred during processing: {str(e)}")

if __name__ == '__main__':
    main()
