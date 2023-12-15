import os
import concurrent.futures
import gzip
import zipfile
import re
import argparse
from functools import partial
from dataclasses import dataclass
from typing import List, Tuple, Any, Optional
import logging

logging.basicConfig(level=logging.INFO)  # Set up logging

# Define a data class for Sequence
@dataclass
class Sequence:
    id: str
    sequence: str
    quality: str = ''

# Function to generate regular expression for primer matching
def generate_regex(primer: str) -> str:
    # Define nucleotide variants for regular expression
    nucleotide_variants = {
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
        variants = nucleotide_variants.get(base, [base])  # Get the variants for the base
        regex_pattern += '[' + ''.join(variants) + ']' if len(variants) > 1 else base  # Construct the regular expression pattern
    return regex_pattern

# Function to get reverse complement of a DNA sequence
def reverse_complement(dna_sequence: str) -> str:
    # Define mapping for complement bases
    complement_bases = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    reverse_complement_sequence = ''.join(complement_bases.get(base, base) for base in reversed(dna_sequence))  # Get the reverse complement
    return reverse_complement_sequence

# Function to check if a sequence contains both forward and reverse primers
def extract_trimmed_sequence(sequence_obj: Sequence, forward_primer_regex: Any, reverse_primer_regex: Any) -> Optional[Sequence]:
    forward_match = re.search(forward_primer_regex, sequence_obj.sequence)
    reverse_match = re.search(reverse_primer_regex, reverse_complement(sequence_obj.sequence))
    
    if forward_match and reverse_match:
        # Calculate start and end positions for trimming
        start_pos = forward_match.end()
        end_pos = len(sequence_obj.sequence) - reverse_match.end()
        return Sequence(sequence_obj.id, sequence_obj.sequence[start_pos:end_pos])
    return None

# Function to check primers for a batch of DNA sequences
def check_batch(batch_of_sequences: List[Sequence], forward_primer_regex: Any, reverse_primer_regex: Any) -> List[str]:
    # Return the ids of sequences that have both primers
    return [extract_trimmed_sequence(sequence_obj, forward_primer_regex, reverse_primer_regex) for sequence_obj in batch_of_sequences]

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

# Main function to process the input DNA file and return the matching sequence ids
def process_file_with_primers(input_path: str, forward_primer: str, reverse_primer: str, num_threads: int) -> List[Sequence]:
    logging.info(f"Processing file: {input_path}")
    sequences_with_id = []

    # Validate inputs
    if not os.path.isfile(input_path):
        logging.error("The provided input path does not point to a file.")
        raise FileNotFoundError("The provided input path does not point to a file.")

    if not all((isinstance(forward_primer, str), isinstance(reverse_primer, str))):
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

    is_fastq = any(input_path.endswith(extension) for extension in ['.fastq', '.fq', '.fastq.gz', '.fq.gz'])
    current_id = None
    current_seq = ""
    current_quality = ""
    for line in file_lines:
        line = line.strip()
        if line.startswith(">") or line.startswith("@"):
            if current_id and current_seq:
                sequences_with_id.append(Sequence(current_id, current_seq, current_quality))
            current_id = line[1:]
            current_seq = ""
            current_quality = ""
        elif is_fastq and line.startswith("+"):
            current_quality = next(file_lines).strip()
        else:
            current_seq += line

    if current_id and current_seq:
        sequences_with_id.append(Sequence(current_id, current_seq, current_quality))

    logging.info("Generating regular expressions for primers.")
    forward_primer_regex = generate_regex(forward_primer)  # Generate regular expression for forward primer
    reverse_primer_regex = generate_regex(reverse_primer)  # Generate regular expression for reverse primer

    logging.info("Checking sequences for primers.")
    # Calculate the batch size for dividing the sequences between threads
    # This is done using ceiling division to ensure all sequences are covered
    batch_size = -(-len(sequences_with_id) // num_threads)  # Ceiling division

    # Using functools.partial, we create a new function that has the same
    # functionality as check_batch but also has the regular expressions for
    # the primers pre-filled
    check_batch_with_primers = partial(check_batch, forward_primer_regex=forward_primer_regex, reverse_primer_regex=reverse_primer_regex)  # Create a new function with regular expressions pre-filled

    # Split the list of sequences into batches for processing in parallel
    sequence_batches = [sequences_with_id[i:i+batch_size] for i in range(0, len(sequences_with_id), batch_size)]  # Divide the sequences into batches

    # Process each batch of sequences in parallel using concurrent.futures
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:  # Create a process pool executor
        results = executor.map(check_batch_with_primers, sequence_batches)  # Apply the function to each batch

    # Combine the results from all batches into a single list
    matching_sequences = [seq for batch in results for seq in batch if seq is not None]  # Collect the matching Sequence object

    logging.info("Done processing file.")
    return matching_sequences  # Return the matching Sequence object

def main():
    parser = argparse.ArgumentParser(description='Process DNA sequences with primers.')
    parser.add_argument('-i', '--input_path', type=str, help='Path to the input file', required=True)
    parser.add_argument('-f', '--forward_primer', type=str, help='Forward primer sequence', required=True)
    parser.add_argument('-r', '--reverse_primer', type=str, help='Reverse primer sequence', required=True)
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads (default: 4)')

    args = parser.parse_args()  # Parse the command line arguments

    try:
        matching_sequences = process_file_with_primers(args.input_path, args.forward_primer, args.reverse_primer, args.num_threads)  # Process the file with primers
        for seq in matching_sequences:
            print(f">{seq.id}\n{seq.sequence}")  # Print the matching sequence ids
    except Exception as e:
        logging.error(f"An error occurred during processing: {str(e)}")

if __name__ == '__main__':
    main()
