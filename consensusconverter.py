import argparse
import logging
import os
from functools import partial
from dataclasses import dataclass
from itertools import product, groupby
from typing import List
import tempfile
import concurrent.futures

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    def __hash__(self):
        return hash((self.id, self.sequence))

def delineate(string: str) -> List[str]:
    # A.................Adenine
    # C.................Cytosine
    # G.................Guanine
    # T (or U)..........Thymine (or Uracil)
    # R.................A or G
    # Y.................C or T
    # S.................G or C
    # W.................A or T
    # K.................G or T
    # M.................A or C
    # B.................C or G or T
    # D.................A or G or T
    # H.................A or C or T
    # V.................A or C or G
    # N.................any base
    # . or -............gap
    conversion = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "T": ["T"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "T", "G", "C"],
    }
    res = []
    for sub in [zip(conversion.keys(), char) for char in product(*conversion.values())]:
        tmp = string
        for repls in sub:
            tmp = tmp.replace(*repls)
        res.append(tmp)
    return list(set(res))

def read_sequences_from_fasta(file_path: str) -> List[Seq]:
    sequences = []
    with open(file_path, 'rt') as fin:
        faiter = (x[1] for x in groupby(fin, lambda line: line[0] == ">"))
        for header in faiter:
            headerstr = next(header).strip()
            name = headerstr[1:]  # Removing ">" character
            seq = "".join(s.strip().upper() for s in next(faiter))
            sequences.append(Seq(name, seq))
    return sequences

def write_sequences_to_fasta(sequences: List[Seq], file_path: str) -> None:
    with open(file_path, 'w') as f:
        for seq in sequences:
            f.write(f">{seq.id}\n{seq.sequence}\n")

def delinate_chunk(sequences: List[Seq]) -> List[Seq]:
    delinated_sequences = []
    for seq in sequences:
        delinated_seqs = delineate(seq.sequence)
        if len(delinated_seqs) == 1:
            delinated_sequences.append(Seq(seq.id, delinated_seqs[0]))
        else:
            max_digits = len(str(len(delinated_seqs)))
            delinated_sequences.extend(
                Seq(f"{seq.id}_{str(index).zfill(max_digits)}", delinated_seq)
                for index, delinated_seq in enumerate(delinated_seqs, start=1)
            )
    return delinated_sequences

def delinate_fasta(input_file: str, output_file: str, num_threads: int) -> None:
    logging.info(f"Reading sequences from {input_file}")
    with tempfile.NamedTemporaryFile(delete=False, mode='w+t', dir='.') as temp_file:
        input_sequences = read_sequences_from_fasta(input_file)
        write_sequences_to_fasta(input_sequences, temp_file.name)

        total_sequences = len(input_sequences)
        chunk_size = max(1, total_sequences // num_threads)
        seq_chunks = [input_sequences[i:i + chunk_size] for i in range(0, total_sequences, chunk_size)]
        logging.info(f"Delinating sequences using {num_threads} threads")
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
            func = partial(delinate_chunk)
            futures = [executor.submit(func, chunk) for chunk in seq_chunks]
            delinated_sequences = []
            for future in concurrent.futures.as_completed(futures):
                delinated_sequences.extend(future.result())
    logging.info(f"Writing delinated sequences to {output_file}")
    write_sequences_to_fasta(delinated_sequences, output_file)

def main():
    parser = argparse.ArgumentParser(description='Delinating sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')

    args = parser.parse_args()

    if not any(ext in args.input_file for ext in ['.fasta', '.fa', '.fna']):
        raise ValueError(f"Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna).")

    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Adjusting thread count to be between 1 and {max_threads}.")
        args.num_threads = min(max(args.num_threads, 1), max_threads)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    delinate_fasta(args.input_file, args.output_file, args.num_threads)

if __name__ == "__main__":
    main()
