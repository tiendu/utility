import os
import re
import gzip
import logging
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial

logging.basicConfig(level=logging.INFO)

def open_file(filename, mode):
    """Helper function to open files based on their extension."""
    try:
        if filename.endswith(".gz"):
            return gzip.open(filename, mode + "t")
        else:
            return open(filename, mode)
    except gzip.BadGzipFile:
        raise Exception(f"Error opening gzip file {filename}. It might be corrupted or not a valid gzip file.")
    except Exception as e:
        raise Exception(f"Error opening file {filename}. Error: {str(e)}")

def interleave_files(prefix, grouped_files, path, pattern):
    sorted_files = sorted(grouped_files, key=lambda f: re.search(pattern, f).group())
    r1_file, r2_file = sorted_files
    output_file = os.path.join(path, f"{prefix}_combined.fastq.gz")

    with open_file(os.path.join(path, r1_file), 'r') as r1, open_file(os.path.join(path, r2_file), 'r') as r2, gzip.open(output_file, 'wt') as out:
        while True:
            r1_lines = [r1.readline() for _ in range(4)]
            r2_lines = [r2.readline() for _ in range(4)]

            # Check if we've reached the end of either file
            if not r1_lines[0] or not r2_lines[0]:
                break

            # Write interleaved reads to output
            out.writelines(r1_lines)
            out.writelines(r2_lines)

    logging.info(f"Interleaved {r1_file} and {r2_file} into {output_file.split('/')[-1]}")

def deinterleave_files(prefix, file, path, pattern):
    input_file = os.path.join(path, file)
    
    # Extract the character inside the square brackets from the pattern
    suffix = re.search(r'\[(.)\]', pattern).group(1)

    r1_output_file = os.path.join(path, f"{prefix}_{suffix}1.fastq.gz")
    r2_output_file = os.path.join(path, f"{prefix}_{suffix}2.fastq.gz")

    with open_file(input_file, 'r') as inp, gzip.open(r1_output_file, 'wt') as r1_out, gzip.open(r2_output_file, 'wt') as r2_out:
        while True:
            r1_lines = [inp.readline() for _ in range(4)]
            r2_lines = [inp.readline() for _ in range(4)]
            
            if not r1_lines[0] or not r2_lines[0]:
                break

            r1_out.writelines(r1_lines)
            r2_out.writelines(r2_lines)

    logging.info(f"Deinterleaved {file} into {r1_output_file.split('/')[-1]} and {r2_output_file.split('/')[-1]}")

def is_interleaved(file_path):
    try:
        with open_file(file_path, 'r') as f:
            # Check first 5 reads
            for _ in range(5):
                r1_header = f.readline().split()[0]
                f.readline()  # skip sequence
                f.readline()  # skip + sign
                f.readline()  # skip quality scores

                r2_header = f.readline().split()[0]
                f.readline()  # skip sequence
                f.readline()  # skip + sign
                f.readline()  # skip quality scores

                # Compare headers
                if r1_header != r2_header:
                    return False
        return True
    except Exception as e:
        logging.error(f"Error reading file {file_path}. Error: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Interleave or deinterleave FASTQ files.")
    parser.add_argument("-p", "--path", required=True, help="Directory path containing the FASTQ files.")
    parser.add_argument("-r", "--pattern", default="_R[12]", help="Regex pattern to match the FASTQ files. Default is _R[12].")
    parser.add_argument("-n", "--num_threads", type=int, default=4, help="Number of threads to use. Default is 4.")
    parser.add_argument("-m", "--mode", choices=["interleave", "deinterleave"], required=True, help="Choose between interleaving or deinterleaving the FASTQ files.")
    
    args = parser.parse_args()

    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Adjusting thread count to be between 1 and {max_threads}.")
        args.num_threads = min(max(args.num_threads, 1), max_threads)

    # List all files in the provided directory
    files = os.listdir(args.path)

    # Group the files by the prefix before the pattern
    file_groups = {}
    for filename in files:
        match = re.search(args.pattern, filename)
        if match:
            prefix = filename[:match.start()]
            if prefix not in file_groups:
                file_groups[prefix] = []
            file_groups[prefix].append(filename)

    with ProcessPoolExecutor(max_workers=args.num_threads) as executor:
        if args.mode == "interleave":
            for prefix, grouped_files in file_groups.items():
                if len(grouped_files) != 2:
                    logging.warning(f"{prefix} does not have both R1 and R2 files. Skipping...")
                    continue
                executor.submit(interleave_files, prefix, grouped_files, args.path, args.pattern)
        elif args.mode == "deinterleave":
            combined_files = [f for f in files if not re.search(args.pattern, f) and is_interleaved(os.path.join(args.path, f))]
            for file in combined_files:
                prefix = file.split("_combined")[0]
                worker = partial(deinterleave_files, prefix=prefix, file=file, path=args.path, pattern=args.pattern)
                executor.submit(worker)

if __name__ == "__main__":
    main()
