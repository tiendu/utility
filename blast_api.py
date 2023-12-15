import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import logging
import os
import random
import time
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from functools import partial
from itertools import groupby
from multiprocessing import Manager, Lock
from typing import List, Dict, Tuple
from urllib.parse import quote
import tempfile
import requests


@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __hash__(self):
        return hash((self.id, self.sequence, self.quality))

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

def hash_sequence(sequence: str) -> str:
    """Returns a SHA-256 hash of a sequence."""
    return hashlib.sha256(sequence.encode()).hexdigest()

def deduplicate_chunk(sequences: List[Seq], seen_hashes: Dict[str, bool], lock: Lock) -> List[Seq]:
    logging.info(f"Processing a chunk with {len(sequences)} sequences.")
    sequences.sort(key=lambda s: len(s.sequence), reverse=True)  # Sort by length
    unique_seqs = []

    for current_seq in sequences:
        seq_hash = hash_sequence(current_seq.sequence)

        # Check if sequence is already seen using its hash
        if seq_hash in seen_hashes:
            continue

        # Check if sequence is contained in any of the unique sequences
        is_contained = any(current_seq.sequence in unique_seq.sequence for unique_seq in unique_seqs)
        if is_contained:
            continue

        unique_seqs.append(current_seq)
        with lock:
            seen_hashes[seq_hash] = True

    logging.info(f"Deduplicated chunk to {len(unique_seqs)} unique sequences.")
    return unique_seqs

def recursive_deduplication(input_file: str, file_type: str, num_threads: int, seen_hashes_shared: Manager().dict, lock: Lock) -> List[Seq]:
    deduped_seqs_all = set()

    while True:
        sequences = read_sequences_from_file(input_file, file_type)

        if not sequences:
            break

        logging.info(f"Initial number of sequences: {len(sequences)}")

        # Split sequences into chunks.
        total_sequences = len(sequences)
        chunk_size = max(1, total_sequences // num_threads)
        seq_chunks = [sequences[i:i + chunk_size] for i in range(0, total_sequences, chunk_size)]

        deduped_seqs = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            func = partial(deduplicate_chunk, seen_hashes=seen_hashes_shared, lock=lock)
            futures = [executor.submit(func, chunk) for chunk in seq_chunks if chunk]

            for future in concurrent.futures.as_completed(futures):
                deduped_seqs.extend(future.result())

        # Add deduplicated sequences to the main set
        deduped_seqs_all.update(set(deduped_seqs))

        # If no sequences were deduplicated in this iteration, we're done
        if len(deduped_seqs_all) == len(sequences):
            break

        # Otherwise, shuffle the deduplicated sequences and write them to the input temp file
        sequences = list(deduped_seqs_all)
        random.shuffle(sequences)
        write_sequences_to_file(sequences, input_file)

    return list(deduped_seqs_all)

def cleanup_temp_file(file_path: str) -> None:
    try:
        os.remove(file_path)
    except Exception as e:
        logging.warning(f"Unable to delete temporary file {file_path}. Error: {e}")

def deduplicate_fasta(input_file: str, file_type: str, output_file: str, num_threads: int) -> None:
    with tempfile.NamedTemporaryFile(delete=False, mode='w+t', dir='.') as temp_file:
        input_sequences = read_sequences_from_file(input_file, file_type)
        write_sequences_to_file(input_sequences, temp_file.name)

        with Manager() as manager:
            seen_hashes_shared = manager.dict()
            lock = manager.Lock()
            deduped_seqs = recursive_deduplication(temp_file.name, file_type, num_threads, seen_hashes_shared, lock)
        write_sequences_to_file(deduped_seqs, output_file)
        logging.info(f"Wrote {len(deduped_seqs)} unique sequences to {output_file}.")

        cleanup_temp_file(temp_file.name)

MAX_RETRIES = 10
RATE_LIMIT = 1

def poll_blast_results(rid: str) -> bool:
    retries = 0
    sleep_time = RATE_LIMIT

    while retries < MAX_RETRIES:
        try:
            url = f"https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}"
            response = requests.get(url)
            content = response.text

            if "Status=WAITING" in content:
                retries += 1
                time.sleep(sleep_time)
                sleep_time *= 2  # Exponential backoff
                continue
            if "Status=FAILED" in content:
                raise ValueError("Search failed; please report to blast-help@ncbi.nlm.nih.gov.")
            if "Status=UNKNOWN" in content:
                raise ValueError("Search expired.")
            if "Status=READY" in content:
                return "ThereAreHits=yes" in content

        except requests.exceptions.ConnectionError as e:
            logging.warning(f"Connection error encountered: {e}. Retrying...")
            retries += 1
            time.sleep(sleep_time)
            sleep_time *= 2  # Exponential backoff
            continue

    logging.error(f"Max retries reached for RID {rid}. Moving to the next sequence.")
    return False

def initiate_blast_search(sequence: str, database: str, service: str) -> Tuple[str, int]:
    url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    query = quote(sequence)
    data = {
        "CMD": "Put",
        "PROGRAM": service,
        "DATABASE": database,
        "QUERY": query,
        "FORMAT_TYPE": "XML",  # Request XML format
        "HITLIST_SIZE": "5"   # Limit results to 5 hits
    }

    try:
        response = requests.post(url, data=data, timeout=30)
        response.raise_for_status()
        content = response.text

        rid_match = re.search(r"RID = (\S+)", content)
        rtoe_match = re.search(r"RTOE = (\d+)", content)
        
        if not rid_match or not rtoe_match:
            raise ValueError("Failed to initiate BLAST search")
        
        rid = rid_match.group(1)
        rtoe = int(rtoe_match.group(1))
        
        return rid, rtoe
    except requests.RequestException as e:
        raise ValueError(f"HTTP request failed: {e}")

def retrieve_blast_results(rid: str) -> str:
    url = f"https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"
    response = requests.get(url)
    return response.text

def perform_blast(seq_obj: Seq, database: str, service: str) -> str:
    rid, rtoe = initiate_blast_search(seq_obj.sequence, database, service)

    time.sleep(rtoe)

    if poll_blast_results(rid):
        results = retrieve_blast_results(rid)
        return results
    else:
        logging.warning(f"No hits found for sequence {seq_obj.id}")
        return None

def main():
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences and perform BLAST.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-d', '--database', default='nr', help='BLAST database to use.')
    parser.add_argument('-s', '--service', default='blastn', help='BLAST service to use.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')

    args = parser.parse_args()

    file_type = ''
    if any(ext in args.input_file for ext in ['.fastq', '.fq']):
        file_type = 'FASTQ'
    elif any(ext in args.input_file for ext in ['.fasta', '.fa', '.fna']):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).')

    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Adjusting thread count to be between 1 and {max_threads}.")
        args.num_threads = min(max(args.num_threads, 1), max_threads)

    dedup_filename = 'deduped_' + os.path.basename(args.input_file)
    deduplicate_fasta(args.input_file, file_type, dedup_filename, args.num_threads)
    deduped_seqs = read_sequences_from_file(dedup_filename, file_type)

    # Open TSV file for writing results
    with open(args.output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        # Write header to TSV
        writer.writerow(['Hit ID', 'Hit Definition', 'Score', 'E-Value', 'Percent Identity'])
        # BLAST deduplicated sequences and write results to TSV
        results = set()
        for seq in deduped_seqs:
            try:
                result = perform_blast(seq, args.database, args.service)
            except ValueError as e:
                logging.warning(f"Error encountered for sequence {seq.id}: {e}")
            if result:
                root = ET.fromstring(result)
                for hit in root.findall(".//Hit"):
                    hit_id = hit.find("Hit_id").text
                    hit_def = hit.find("Hit_def").text
                    hit_score = hit.find(".//Hsp_score").text
                    hit_evalue = hit.find(".//Hsp_evalue").text
                    identity = int(hit.find(".//Hsp_identity").text)
                    align_len = int(hit.find(".//Hsp_align-len").text)
                    percent_identity = (identity / align_len) * 100
                    
                    formatted_result = f"{seq.id}\t{hit_id}\t{hit_def}\t{hit_score}\t{hit_evalue}\t{percent_identity:.2f}%"
                    print(formatted_result)
                    if formatted_result not in results:
                        writer.writerow(formatted_result.split('\t'))
                        results.add(formatted_result)

if __name__ == "__main__":
    main()
