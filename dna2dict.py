from html import unescape
import sys
import struct
import re

def remove_html_tags(text):
    '''Remove HTML/XML tags from the text.'''
    clean_text = re.sub(r'<[^>]*>', '', text)

    return clean_text

def clean_notes(text):
    # Regular expression to capture tag content
    pattern = re.compile(r'<(\w+)>(.*?)</\1>', re.DOTALL)

    notes = []
    
    # Extract and decode content
    for match in pattern.finditer(text):
        tag, content = match.groups()
        # Decode HTML entities
        decoded_content = unescape(content)
        decoded_content = re.sub(r'<.*?>', '', decoded_content)
        notes.append(f'"{tag}": "{decoded_content}"')

    return '; '.join(notes)

def unpack_data(file_obj, size, mode):
    '''Unpack data from the file object.'''
    return struct.unpack('>' + mode, file_obj.read(size))[0]

def parse_snapgene_file(filepath=None, file_obj=None):
    '''Return a dictionary containing the data from a .dna file created with SnapGene.'''
    if filepath is not None:
        file_obj = open(filepath, 'rb')

    if file_obj.read(1) != b'\t':
        raise ValueError('Invalid SnapGene file format!')

    # Read the document properties
    file_length = unpack_data(file_obj, 4, 'I')
    file_title = file_obj.read(8).decode('ascii')

    if file_length != 14 or file_title != 'SnapGene':
        raise ValueError('Invalid SnapGene file format!')

    snapgene_data = {
        'isDNA': unpack_data(file_obj, 2, 'H'),
        'exportVersion': unpack_data(file_obj, 2, 'H'),
        'importVersion': unpack_data(file_obj, 2, 'H'),
        'features': []
    }
    strand_mapping = {'0': '.', '1': '+', '2': '-', '3': '='}

    while True:
        block_type_byte = file_obj.read(1)

        if block_type_byte == b'':
            # End of file
            break

        block_size = unpack_data(file_obj, 4, 'I')

        if ord(block_type_byte) == 0:
            # Read the DNA sequence and its properties
            sequence_properties = unpack_data(file_obj, 1, 'b')
            snapgene_data['dna'] = {
                'topology': 'circular' if sequence_properties & 0x01 else 'linear',
                'strandedness': 'double' if sequence_properties & 0x02 > 0 else 'single',
                'damMethylated': sequence_properties & 0x04 > 0,
                'dcmMethylated': sequence_properties & 0x08 > 0,
                'ecoKIMethylated': sequence_properties & 0x10 > 0,
                'length': block_size - 1
            }
            snapgene_data['sequence'] = file_obj.read(block_size - 1).decode('ascii')
        elif ord(block_type_byte) == 6:
            # Read the notes
            block_content = file_obj.read(block_size).decode('utf-8')
            snapgene_data['notes'] = parse_notes(block_content)
        elif ord(block_type_byte) == 10:
            # Read the features
            features_content = file_obj.read(block_size).decode('utf-8')
            snapgene_data['features'].extend(parse_features(features_content, strand_mapping))
        else:
            # Ignore the block
            file_obj.read(block_size)

    if filepath is not None:
        file_obj.close()

    return snapgene_data

def parse_notes(block_content):
    '''Parse notes from the block content.'''
    notes = {}
    # Simple regex to extract key-value pairs in the notes section
    note_pattern = re.compile(r'<(\w+)>(.*?)</\1>', re.DOTALL)

    for match in note_pattern.findall(block_content):
        key, value = match
        notes[key] = clean_notes(value)

    return notes

def parse_features(features_content, strand_mapping):
    '''Parse features from the features content.'''
    features = []
    feature_pattern = re.compile(r'<Feature(.*?)>(.*?)</Feature>', re.DOTALL)
    segment_pattern = re.compile(r'<Segment(.*?)>', re.DOTALL)
    qualifier_pattern = re.compile(r'<Q name="([^"]+)">(.*?)</Q>', re.DOTALL)

    for feature_match in feature_pattern.findall(features_content):
        feature_attributes = feature_match[0].strip()
        feature_body = feature_match[1].strip()

        # Extract attributes like type, name, and directionality
        feature_type = re.search(r'type="([^"]+)"', feature_attributes).group(1) if re.search(r'type="([^"]+)"', feature_attributes) else None
        feature_name = re.search(r'name="([^"]+)"', feature_attributes).group(1) if re.search(r'name="([^"]+)"', feature_attributes) else None
        directionality = re.search(r'directionality="([^"]+)"', feature_attributes)
        strand = strand_mapping[directionality.group(1)] if directionality else strand_mapping['0']
        segments = []

        for segment_match in segment_pattern.findall(feature_body):
            segment_range = re.search(r'range="([^"]+)"', segment_match).group(1)
            start, end = [int(x) for x in segment_range.split('-')]
            color = re.search(r'color="([^"]+)"', segment_match).group(1)
            segments.append({'start': start, 'end': end, 'color': color})

        qualifiers = {}

        for qualifier_match in qualifier_pattern.findall(feature_body):
            qualifier_name, qualifier_value = qualifier_match
            qualifiers[qualifier_name] = remove_html_tags(qualifier_value.strip())

        if 'label' not in qualifiers:
            qualifiers['label'] = feature_name

        if 'note' not in qualifiers:
            qualifiers['note'] = []

        if not isinstance(qualifiers['note'], list):
            qualifiers['note'] = [qualifiers['note']]

        qualifiers['note'].append('color: ' + segments[0]['color'])
        features.append({
            'start': min(segment['start'] - 1 for segment in segments),
            'end': max(segment['end'] for segment in segments),
            'strand': strand,
            'type': feature_type,
            'name': feature_name,
            'color': segments[0]['color'],
            'textColor': 'black',
            'segments': segments,
            'row': 0,
            'isOrf': False,
            'qualifiers': qualifiers
        })

    return features

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python extract_features.py <path_to_dna_file>')
        sys.exit(1)

    dna_file_path = sys.argv[1]
    result = parse_snapgene_file(filepath=dna_file_path)
    print(result)
