def dna_to_bits(dna: str) -> list:
    nu_to_bits = {
        'A': 0b0001, 'T': 0b0010, 'G': 0b0100, 'C': 0b1000,
        'W': 0b0011, 'R': 0b0101, 'Y': 0b1010, 'S': 0b1100,
        'K': 0b0110, 'M': 0b1001, 'B': 0b1110, 'D': 0b0111,
        'H': 0b1011, 'V': 0b1101, 'N': 0b1111
    }
    return [nu_to_bits[nu.upper()] for nu in dna]

def compute_lps(pattern: list) -> list:
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

def kmp_search(text: str, pattern: str, modifier=None, tolerance=0) -> list:
    def match(a, b):
        if isinstance(a, int) and isinstance(b, int):
            return (a & b) != 0
        return a == b

    if modifier:
        text = modifier(text)  # Encode text if modifier is provided
        pattern = modifier(pattern)  # Encode pattern if modifier is provided
    else:
        text = list(text)  # Treat as a list of characters for direct comparison
        pattern = list(pattern)

    if tolerance >= len(pattern):
        raise ValueError("Tolerance cannot be greater than or equal to the length of the pattern.")

    lps = compute_lps(pattern)  # Compute LPS array
    occurrences = []  # List to store occurrences of the pattern

    for i in range(len(text)):  # Iterate over each character in text
        mismatch = 0  # Initialize mismatch count for each starting position
        j = 0  # Reset pattern index for each starting position

        # Continue matching as long as we're within bounds
        while j < len(pattern) and (i + j) < len(text):
            if match(text[i + j], pattern[j]):
                # Move to the next character if it's a match
                j += 1
            else:
                # Count the mismatch
                mismatch += 1
                if mismatch > tolerance:
                    break  # Exit if mismatches exceed tolerance
                j += 1

        # If j reached the end of the pattern, we found a match
        if j == len(pattern) and mismatch <= tolerance:
            occurrences.append(i)  # Append the starting index of the match

    return occurrences

# Example usage
text = "ATGCGATACGCTAGCTAGCATCGAT"
pattern = "ATACN"
result = kmp_search(text, pattern, dna_to_bits, tolerance=int(len(pattern)*0.2))

print(f"Pattern found at indices: {result}")
