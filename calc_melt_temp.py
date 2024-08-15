def calculate_melting_temperature(sequence: str, salt_conc: float=50, dna_conc: float=0.5, is_rna: bool=False) -> float:
    '''
    Calculate the melting temperature of a DNA or RNA sequence using the nearest-neighbor thermodynamic model.
    Source: https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/documents/367/000/meltingtemp1.pdf

    Parameters:
    - sequence (str): DNA or RNA sequence.
    - salt_conc (float): Salt concentration in millimolar (mM). Default is 50 mM.
    - dna_conc (float): DNA or RNA concentration in micromolar (µM). Default is 0.5 µM.
    - is_rna (bool): True if the sequence is RNA, False if it is DNA.

    Returns:
    - float: Melting temperature in degrees Celsius.
    '''
    # Nearest-neighbor parameters for DNA
    nn_params_dna = {
        ('A', 'A'): (-9.1, -24.0),
        ('A', 'T'): (-8.6, -23.9),
        ('T', 'A'): (-6.0, -16.9),
        ('C', 'A'): (-5.8, -12.9),
        ('G', 'T'): (-6.5, -17.3),
        ('C', 'T'): (-7.8, -20.8),
        ('G', 'A'): (-5.6, -13.5),
        ('C', 'G'): (-11.9, -27.8),
        ('G', 'C'): (-11.1, -26.7),
        ('G', 'G'): (-11.0, -26.6)
    }

    # Nearest-neighbor parameters for RNA
    nn_params_rna = {
        ('T', 'T'): (-6.6, -18.4),
        ('T', 'A'): (-5.7, -15.5),
        ('A', 'T'): (-8.1, -22.6),
        ('G', 'T'): (-10.5, -27.8),
        ('C', 'A'): (-10.2, -26.2),
        ('G', 'A'): (-7.6, -19.2),
        ('C', 'T'): (-13.3, -35.5),
        ('G', 'C'): (-8.0, -19.4),
        ('C', 'G'): (-14.2, -34.9),
        ('C', 'C'): (-12.2, -29.7)
    }

    # Choose parameters based on sequence type
    nn_params = nn_params_rna if is_rna else nn_params_dna

    # Initialize parameters
    delta_H = 0.0
    delta_S = 0.0
    length = len(sequence)

    # Calculate the sum of enthalpy and entropy changes
    for i in range(length - 1):
        pair = (sequence[i], sequence[i + 1])
        if pair in nn_params:
            delta_H += nn_params[pair][0]
            delta_S += nn_params[pair][1]

    # Add initiation parameters
    initiation_H = 0.0
    initiation_S = -10.8
    delta_H += initiation_H
    delta_S += initiation_S

    # Constants
    R = 1.99  # gas constant in cal/(mol*K)
    C_molar = dna_conc * 1e-6  # Convert µM to M
    Na_conc = salt_conc * 1e-3  # Convert mM to M

    # Melting temperature calculation
    Tm = (1000 * delta_H) / (delta_S + R * math.log(C_molar / 4)) - 273.15 + 16.6 * math.log10(Na_conc)

    return Tm

def check_hairpin(primer: str, min_stem_length: int=4, max_loop_length: int=10) -> bool:
    '''
    Check for hairpin formation in a primer sequence.

    Parameters:
    - primer (str): Primer sequence.
    - min_stem_length (int): Minimum length of the stem for a hairpin. Default is 4.
    - max_loop_length (int): Maximum length of the loop allowed. Default is 10.

    Returns:
    - bool: True if hairpin is detected, False otherwise.
    '''
    primer_length = len(primer)

    def is_self_complementary(seq1: str, seq2: str) -> bool:
        '''Check if two sequences are self-complementary.'''
        complement = str.maketrans('ATGC', 'TACG')
        return seq1.translate(complement) == seq2

    # Check for hairpin formation
    for stem_length in range(min_stem_length, primer_length // 2):
        for start in range(primer_length - 2 * stem_length):
            stem1 = primer[start:start + stem_length]
            stem2 = primer[start + stem_length:start + 2 * stem_length]
            loop_length = start + 2 * stem_length - primer_length

            # Check if the ends are self-complementary
            if loop_length <= max_loop_length and is_self_complementary(stem1, stem2[::-1]):
                return True

    return False

def check_uniqueness(primer: str, target_sequence: str) -> bool:
    '''
    Check if a primer is unique within the target sequence using regex.

    Parameters:
    - primer (str): Primer sequence.
    - target_sequence (str): Target DNA sequence.

    Returns:
    - bool: True if primer is unique (appears only once), False otherwise.
    '''
    # Create a regex pattern to find all occurrences of the primer
    pattern = re.compile(re.escape(primer))

    # Find all matches in the target sequence
    matches = list(pattern.finditer(target_sequence))

    # Check if primer is found more than once
    return len(matches) == 1

def generate_flanking_primers(sequence: str, target_start: int, target_end: int, primer_length: int=20, tm_min: float=50, tm_max: float=60) -> None:
    '''
    Generate primers that flank a specific target region in a given DNA sequence, checking for Tm, uniqueness, and hairpin formation.

    Parameters:
    - sequence (str): DNA sequence.
    - target_start (int): Start position of the target region (1-based index).
    - target_end (int): End position of the target region (1-based index).
    - primer_length (int): Length of the primers to be generated. Default is 20.
    - tm_min (float): Minimum melting temperature for primers. Default is 50°C.
    - tm_max (float): Maximum melting temperature for primers. Default is 60°C.

    Returns:
    - list: List of potential primers meeting the criteria.
    '''
    primers = []

    # Generate upstream primers
    for i in range(target_start - primer_length):
        primer = sequence[i:i + primer_length]
        tm = calculate_melting_temperature(primer)
        if tm_min <= tm <= tm_max and not check_hairpin(primer):
            if check_uniqueness(primer, sequence):
                primers.append({
                    'sequence': primer,
                    'position': 'upstream',
                    'start': i,
                    'end': i + primer_length - 1,
                    'Tm': tm
                })

    # Generate downstream primers
    for i in range(target_end, len(sequence) - primer_length + 1):
        primer = sequence[i:i + primer_length]
        tm = calculate_melting_temperature(primer)
        if tm_min <= tm <= tm_max and not check_hairpin(primer):
            if check_uniqueness(primer, sequence):
                primers.append({
                    'sequence': primer,
                    'position': 'downstream',
                    'start': i,
                    'end': i + primer_length - 1,
                    'Tm': tm
                })

    for primer in primers:
        print(f'Primer: {primer["sequence"]}, Position: {primer["position"]}, Start: {primer["start"]}, End: {primer["end"]}, Tm: {primer["Tm"]:.2f} °C')
