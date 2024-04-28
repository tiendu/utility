down_arrow = "↓"
def print_changes(before, after, width=5):
    # Print the sequences and annotations
    for i in range(0, len(before), width):
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
        print(f'{i+1} {after_chunk}')

# Example usage
before_seq = "ATGCACCTG"
after_seq = "ATGCTCCTG"
print_changes(before_seq, after_seq)

"""
before_seq = "ATGCACCTG"
after_seq = "ATGCTCCTG"
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
