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

@dataclass()
class Seq:
    id: str
    sequence: str
    def __hash__(self):
        return hash((self.id, self.sequence))

def reverse_translate(sequence: str) -> List[str]:
    aa_vars = {
        'A': ['GCN'],
        'C': ['TGY'],
        'D': ['GAY'],
        'E': ['GAR'],
        'F': ['TTY'],
        'G': ['GGN'],
        'H': ['CAY'],
        'I': ['ATH'],
        'K': ['AAR'],
        'L': ['YTR', 'CTN'],
        'M': ['ATG'],
        'N': ['AAY'],
        'P': ['CCN'],
        'Q': ['CAR'],
        'R': ['CGN', 'MGR'],
        'S': ['TCN', 'AGY'],
        'T': ['ACN'],
        'V': ['GTN'],
        'W': ['TGG'],
        'Y': ['TAY']
    }
    var_list = [aa_vars.get(aa, [aa]) for aa in sequence]
    var_gen = [''.join(var) for var in product(*var_list)]
    return var_gen

def generate_variants(sequence: str) -> List[str]:
    nu_vars = {
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
    var_list = [nu_vars.get(nu, [nu]) for nu in sequence]
    var_gen = [''.join(var) for var in product(*var_list)]
    return var_gen

def read_sequences_from_file(file_path: str) -> List[Seq]:
    sequences = []

    # Determine the appropriate file opening mode based on the file extension.
    opener = gzip.open if file_path.endswith('.gz') else open

    with opener(file_path, 'rt') as fin:
        faiter = (x[1] for x in groupby(fin, lambda line: line[0] == ">"))
        for header in faiter:
            headerstr = next(header).strip()
            name = headerstr[1:]  # Removing ">" character
            seq = "".join(s.strip().upper() for s in next(faiter))
            if "B" in seq.upper() or not seq.isalpha():
                raise ValueError(f"{name} is not comprised of amino acids!")
            sequences.append(Seq(name, seq.upper()))
    return sequences

def write_sequences_to_fasta(sequences: List[Seq], file_path: str) -> None:
    with open(file_path, 'w') as f:
        for seq in sequences:
            f.write(f">{seq.id}\n{seq.sequence}\n")

def delineate_chunk(sequences: List[Seq]) -> List[Seq]:
    delineated_sequences = []
    for seq in sequences:
        reverse_translated = reverse_translate(seq.sequence)
        variants = [generate_variants(i) for i in reverse_translated]
        variants = [j for i in variants for j in i]
        if len(reverse_translated) == 1:
            delineated_sequences.append(Seq(seq.id, variants[0]))
        else:
            max_digits = len(str(len(variants)))
            delineated_sequences.extend(
                Seq(f"{seq.id}_{str(index).zfill(max_digits)}", variant)
                for index, variant in enumerate(variants, start=1)
            )
    return delineated_sequences

def process_file(input_file: str, output_file: str, num_threads: int) -> None:
    logging.info(f"Processing file: {input_file}")
    try:
        input_sequences = read_sequences_from_file(input_file)
    except Exception as e:
        raise ValueError(f'An error occurred while reading the file: {str(e)}')

    # Sort sequences based on length in descending order
    input_sequences.sort(key=lambda x: len(x.sequence), reverse=True)

    total_sequences = len(input_sequences)
    seq_chunks = [[] for _ in range(num_threads)]

    # Distribute sequences in a round-robin fashion
    for i, seq in enumerate(input_sequences):
        seq_chunks[i % num_threads].append(seq)

    logging.info(f"Processing sequences using {num_threads} threads")
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(delineate_chunk, chunk) for chunk in seq_chunks]
        delineated_sequences = []
        for future in concurrent.futures.as_completed(futures):
            delineated_sequences.extend(future.result())
    logging.info(f"Writing delineated sequences to {output_file}")
    write_sequences_to_fasta(delineated_sequences, output_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')

    args = parser.parse_args()

    if not any(ext in args.input_file for ext in ['.fasta', '.fa', '.faa']):
        raise ValueError(f"Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .faa).")

    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Adjusting thread count to be between 1 and {max_threads}.")
        args.num_threads = min(max(args.num_threads, 1), max_threads)

    process_file(args.input_file, args.output_file, args.num_threads)

if __name__ == "__main__":
    main()
