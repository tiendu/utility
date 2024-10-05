
def dna_to_bits(dna: str) -> list:
    """Convert DNA sequence to a list of bitwise representations."""
    nu_to_bits = {
        'A': 0b0001, 'T': 0b0010, 'G': 0b0100, 'C': 0b1000,
        'W': 0b0011, 'R': 0b0101, 'Y': 0b1010, 'S': 0b1100,
        'K': 0b0110, 'M': 0b1001, 'B': 0b1110, 'D': 0b0111,
        'H': 0b1011, 'V': 0b1101, 'N': 0b1111
    }
    return [nu_to_bits[nu.upper()] for nu in dna]

def compute_lps(pattern: list) -> list:
    """Computes the Longest Prefix Suffix (LPS) array for the encoded pattern."""
    lps = [0] * len(pattern)
    length = 0  # Length of the previous longest prefix suffix
    i = 1

    while i < len(pattern):
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1

    return lps

def kmp_search(text: str, pattern: str, modifier=None) -> list:
    """Searches for occurrences of the pattern in the text using the KMP algorithm with an optional modifier for encoding."""

    def match(a, b):
        """Compares two values for a match, considering bitwise and character comparisons."""
        if isinstance(a, int) and isinstance(b, int):
            return (a & b) != 0
        return a == b

    if modifier:
        text = modifier(text)  # Encode text if modifier is provided
        pattern = modifier(pattern)  # Encode pattern if modifier is provided
    else:
        text = list(text)  # Treat as a list of characters for direct comparison
        pattern = list(pattern)

    lps = compute_lps(pattern)  # Compute LPS array
    i = 0  # Index for text
    j = 0  # Index for pattern
    occurrences = []  # List to store occurrences of the pattern

    while i < len(text):
        if j < len(pattern) and match(text[i], pattern[j]):
            # Increment both indices for a match
            i += 1
            j += 1
        else:
            if j == len(pattern):  # Pattern found
                occurrences.append(i - j)  # Append the starting index of the match
                j = lps[j - 1]  # Use LPS to continue searching for the next match
            elif i < len(text) and (j == 0 or not match(text[i], pattern[j])):  # No match
                j = lps[j - 1] if j != 0 else 0  # Use the LPS array to skip comparisons
                if j == 0:
                    i += 1

    return occurrences

# Example usage
text = "ATGCGATACGCTAGCTAGCATCGAT"
pattern = "ATACM"
result = kmp_search(text, pattern, modifier=dna_to_bits)  # Using the modifier to convert DNA to bits

print(f"Pattern found at indices: {result}")
