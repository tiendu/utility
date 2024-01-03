import csv
import argparse
import tempfile
import os
from typing import List, Tuple, Union


def read_input_files(file_paths: List[str]) -> Tuple[List[str], dict]:
    data = {}
    filenames = []

    for idx, file_path in enumerate(file_paths):
        filenames.append(file_path)
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                if len(row) >= 2:
                    key = row[0]
                    value = int(row[1].strip()) if row[1].strip().isdigit() else 0
                    if key not in data:
                        data[key] = [0] * len(file_paths)
                    data[key][idx] = value

    return filenames, data

def write_temp_file(header: List[str], data: dict) -> str:
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        temp_filename = temp_file.name
        temp_file.write("ID\t" + "\t".join(header) + "\n")
        for key, values in data.items():
            temp_file.write(key + "\t" + "\t".join(map(str, values)) + "\n")

    return temp_filename

def merge_files(file_paths: List[str]) -> str:
    filenames, data = read_input_files(file_paths)
    temp_filename = write_temp_file(filenames, data)

    return temp_filename

def append_files(current_table: str, new_table: str) -> str:
    is_header = True
    header = []
    data = {}
    current_row_length = 0

    with open(current_table, 'r', newline='') as csvfile:
         reader = csv.reader(csvfile, delimiter='\t')
         for row in reader:
             if is_header:
                 header = row
                 is_header = False
                 continue
             if len(row) >= 2:
                current_row_length = len(row) - 1
                key = row[0]
                values = [int(value.strip()) if value.strip().isdigit() else 0 for value in row[1:]]
                if key not in data:
                    data[key] = values

    with open(new_table, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if len(row) >= 2:
                key = row[0]
                val = int(row[1].strip()) if row[1].strip().isdigit() else 0
                if key in data:
                    data[key].append(val)
                else:
                    data[key] = [0] * current_row_length
                    data[key].append(val)

    for _, values in data.items():
        if len(values) == current_row_length:
            values.append(0)
    temp_filename = write_temp_file(header, data)

    return temp_filename

def merge_sequences(data: List[List[Union[str, int]]]) -> Tuple[List[str], dict]:
    result = {}
    header = data[0]

    for row in data[1:]:
        row_id = row[0]
        values = [int(value) if value.strip() else 0 for value in row[1:]]
        if row_id not in result:
            result[row_id] = values
        else:
            for i in range(len(values)):
                result[row_id][i] += values[i]
    merged_result = {}

    for current_id, current_value in result.items():
        found = False
        for other_id, other_value in result.items():
            if current_id != other_id and (current_id in other_id or other_id in current_id):
                merged_id = current_id if len(current_id) >= len(other_id) else other_id
                if merged_id not in merged_result:
                    merged_result[merged_id] = [0] * (len(header) - 1)
                for i in range(len(current_value)):
                    merged_result[merged_id][i] += current_value[i]
                found = True
                break
        if not found:
            merged_result[current_id] = current_value

    return header, merged_result

def write_output_file(output_file: str, data: Tuple[List[str], dict]) -> None:
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(data[0])
        for row_id, values in data[1].items():
            writer.writerow([row_id] + values)

def main():
    parser = argparse.ArgumentParser(description='Merge sequences based on substring relationships.')
    parser.add_argument('-i', '--input', nargs='+', help='Path to the input TSV file(s). Provide at least two input files for merging.')
    parser.add_argument('-o', '--output', help='Path to the output TSV file.')
    parser.add_argument('-a', '--append', action='store_true', default=False, help='Append a new table to the current count table.')

    args = parser.parse_args()

    if not args.input or len(args.input) < 2 or not all(os.path.exists(file) for file in args.input):
        raise ValueError("Please provide at least two existing input files.")

    try:
        if not args.append:
            temp_filename = merge_files(args.input)
        else:
            current_table = args.input[0]
            new_table = args.input[1]
            if not all(os.path.exists(file) for file in [current_table, new_table]):
                raise ValueError("Both input files must exist for appending.")
            temp_filename = append_files(current_table, new_table)
        with open(temp_filename, 'r') as temp_file:
            reader = csv.reader(temp_file, delimiter='\t')
            data = [row for row in reader]
        merged_data = merge_sequences(data)
        write_output_file(args.output, merged_data)
    except Exception as e:
        print(f"An error occurred: {str(e)}")

    finally:
        if 'temp_filename' in locals() and os.path.exists(temp_filename):
            os.remove(temp_filename)

if __name__ == "__main__":
    main()
