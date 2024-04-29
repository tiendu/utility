down_arrow = "↓"
def print_changes(before, after, start_pos=1, end_pos=None, width=5):
    # Check if the strings are of equal length
    if len(before) != len(after):
        raise ValueError("Error: strings are not of equal length")
    
    # Set the end position to the length of the strings if not specified
    if end_pos is None:
        end_pos = max(len(before), len(after))

    # Ensure start and end positions are within bounds
    start_pos = max(1, start_pos)
    end_pos = min(end_pos, max(len(before), len(after)))

    # Adjust width if the specified range is within the width of the chunks
    if end_pos - start_pos + 1 < width:
        width = end_pos - start_pos + 1

    # Print the strings and annotations within the specified range
    for i in range(start_pos - 1, end_pos, width):
        before_chunk = before[i:i+width]
        after_chunk = after[i:i+width]
        
        # Construct annotation string for the current chunk
        annotation_chunk = ''
        chunk_length = min(len(before_chunk), len(after_chunk))
        for j in range(chunk_length):
            if before_chunk[j] != after_chunk[j]:
                annotation_chunk += f'{down_arrow}'
            else:
                annotation_chunk += ' '

        # Pad the annotation string if needed to match the width
        if len(annotation_chunk) < width:
            annotation_chunk += ' ' * (width - len(annotation_chunk))

        # Print the chunk with annotations
        print(f'{i+1} {before_chunk}')
        print(f'  {annotation_chunk}')
        print(f'{i+1} {after_chunk}\n')

# Example usage
before = "ATGCACCTGATGCACCTGATGCACCTG"
after = "ATGCTCCTGATGCACCTGATGCACCTG"
print_changes(before, after, width=20)

"""
after = "ATGCACCTG"
before = "ATGCTCCTG"
width = 4
Expected output:
1 ATGC
   ↓    
1 ATGC
5 ACCT
  ↓
5 TCCT
9 G

9 G
"""
