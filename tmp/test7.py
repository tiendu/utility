# Generate k-mers from a sequence
def kmerize(string: str, k: int, modifier=None) -> list:
    if modifier:
        return [modifier(string[i:i+k]) for i in range(len(string) - k + 1)]
    else:
        return [string[i:i+k] for i in range(len(string) - k + 1)]

# Convert a sequence into bit representation
def sequence_to_bits(dna: str) -> list:
    nucl_to_bits = {
        'A': 0b0001,  # A = 1
        'T': 0b0010,  # T = 2
        'G': 0b0100,  # G = 4
        'C': 0b1000,  # C = 8
        'W': 0b0011,  # W = A | T = 3
        'R': 0b0101,  # R = A | G = 5
        'Y': 0b1010,  # Y = C | T = 10
        'S': 0b1100,  # S = G | C = 12
        'K': 0b0110,  # K = G | T = 6
        'M': 0b1001,  # M = A | C = 9
        'B': 0b1110,  # B = C | G | T = 14
        'D': 0b0111,  # D = A | G | T = 7
        'H': 0b1011,  # H = A | C | T = 11
        'V': 0b1101,  # V = A | C | G = 13
        'N': 0b1111   # N = A | C | G | T = 15 (any nucleotide)
    }
    return [nucl_to_bits[nu.upper()] for nu in dna]

# Perform a bitwise comparison of all k-mers from two sequences
def compare_kmer_sets(kmers1: set, kmers2: set) -> int:
    # Compare two k-mers by their bit representation
    def compare_kmers(kmer_bits1: set, kmer_bits2: set):
        return all(bit1 & bit2 for bit1, bit2 in zip(kmer_bits1, kmer_bits2))

    total_matches = 0
    for kmer_bits1 in kmers1:
        # Find the best match in kmers2 for the current kmer_bits1
        best_match = max(compare_kmers(kmer_bits1, kmer_bits2) for kmer_bits2 in kmers2)
        total_matches += best_match
    return total_matches

# Calculate k-mer similarity using bitwise comparison
def calculate_kmer_similarity(seq1: str, seq2: str, k: int) -> float:
    # Generate k-mers from both sequences
    kmers1_bits = kmerize(seq1, k, sequence_to_bits)
    kmers2_bits = kmerize(seq2, k, sequence_to_bits)
    # Perform bitwise comparison and calculate total matches
    total_matches = compare_kmer_sets(kmers1_bits, kmers2_bits)
    # Calculate similarity as the proportion of matching k-mers
    total_kmers = max(len(kmers1_bits), len(kmers2_bits))
    similarity = total_matches / total_kmers if total_kmers > 0 else 0
    return similarity

# Find the best match location for the shorter sequence in the longer sequence
def map_sequences(long_seq: str, 
                  short_seq: str, 
                  similarity_threshold: float=0.0,
                  coverage_threshold: float=0.0,
                  k: int=3) -> list[set[str, float, float]]:
    # Ensure the shorter sequence is the second one
    if len(long_seq) < len(short_seq):
        long_seq, short_seq = short_seq, long_seq
    k = max(len(short_seq) // 5, k)
    coverage = len(short_seq) / len(long_seq)
    results = []
    # Slide the short sequence along the long sequence
    if coverage > coverage_threshold:
        for i in range(len(long_seq) - len(short_seq) + 1):
            # Extract a segment from the long sequence of the same length as the short sequence
            long_segment = long_seq[i:i + len(short_seq)]
            # Calculate similarity between the short sequence and the current segment
            similarity = calculate_kmer_similarity(short_seq, long_segment, k)
            # Update if we find a better match
            if similarity > similarity_threshold:
                similarity = similarity
                position = f'{i + 1}..{i + 1 + len(short_seq)}'
                results.append((position, similarity, coverage))
    return results

# Example sequences
sequence_1 = "tangccatggcta"
sequence_2 = "aygcans"
# Run the mapping and calculate similarity, best position, and coverage
results = map_sequences(sequence_1, sequence_2)
for result in results:
    print(f"Match position: {result[0]}")
    print(f"Similarity score: {result[1]}")
    print(f"Coverage: {result[2]}")
