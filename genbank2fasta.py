import re
import sys

def parse_features(input_file: str) -> dict:
    features = {}
    annotation, id, location = '', None, None
    feature_section = False
    with open(input_file, 'r') as f:
        for line in f:
            line = line.rstrip('\r')
            if 'FEATURES' in line:
                feature_section = True
                continue
            if 'ORIGIN' in line:
                break  
            if feature_section:
                line = line.strip()
                if re.search(r'[0-9]+\.\.[0-9]+', line):
                    id, location = line.split()
                if re.search(r'^/', line):
                    if annotation:
                        features.setdefault(location, {}).setdefault(id, []).append(annotation)
                    annotation = line
                elif annotation and line[0].isalpha():
                    if not re.search(r'[aA-zZ]+ *[0-9]+\.\.[0-9]+', line):
                        annotation += ' ' + line.strip()
                annotation = annotation.replace('/', '')
    return features

def parse_sequence(input_file: str) -> str:
    sequence = ''
    with open(input_file, 'r') as f:
        in_origin = False
        for line in f:
            if 'ORIGIN' in line:
                in_origin = True
                continue
            if in_origin and re.match(r'\s+\d+', line):
                sequence += ''.join(line.strip().split()[1:])
    return sequence.upper()

def extract_sequence(features: dict, sequence: str, output_file: str) -> None:
    reversed_sequence = reverse_complement(sequence)
    with open(output_file, 'w') as f:
        for location, annotations in features.items():
            extract = ''
            header = ''
            if 'complement' in location:
                location = re.search(r'\d+\.\.\d+', location).group()
                start, end = map(int, location.split('..'))
                extract = reversed_sequence[start - 1 : end][::-1]
            elif 'join' in location:
                locations = re.findall(r'\d+\.\.\d+', location)
                for location in locations:
                    start, end = map(int, location.split('..'))
                    extract += sequence[start - 1 : end]
            elif re.match(r'\d+\.\.\d+', location):
                start, end = map(int, location.split('..'))
                extract = sequence[start - 1:end]
            for id, annotation in annotations.items():
                header = id + '|' + location
                header += '|' + '|'.join(annotation)
            f.write(f'>{header}\n{extract}\n')

def reverse_complement(sequence: str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    out = ''
    for i in range(len(sequence), 0, -1):
        out += complement[sequence[i - 1]]
    return out

def main(input_file, output_file):
    features = parse_features(input_file)
    sequence = parse_sequence(input_file)
    extract_sequence(features, sequence, output_file)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: python {sys.argv[0]} input_file output_file')
    else:
        main(sys.argv[1], sys.argv[2])
