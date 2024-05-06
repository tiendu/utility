from dataclasses import dataclass

@dataclass(frozen=True)
class Dna:
    id: str
    sequence: str    
    def __hash__(self):
        return hash((self.id, self.sequence))
    def length(self) -> int:
        return len(self.sequence)
    def gc_content(self) -> float:
        gc_count = sum(1 for base in self.sequence if base in 'GCgc')
        return f'{(gc_count / len(self.sequence)) * 100:.2f}' if self.sequence else 0

def reverse_complement(dna_seq: str) -> str:
    comp_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp_seq = ''.join(comp_bases.get(base, base) for base in reversed(dna_seq))
    return rev_comp_seq

def translate(dna_seq: str) -> str:
    # Dictionary mapping codons to amino acids
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    # Translate DNA sequence into protein sequence
    protein_seq = ''
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        amino_acid = codon_table.get(codon, '')
        if amino_acid == '*':
            break
        protein_seq += amino_acid
    
    return protein_seq

def find_orf(dna_seq: str) -> list:
    start_codons = ['ATG', 'TTG', 'CTG', 'GTG']
    stop_codons = ['TAA', 'TAG', 'TGA']

    orfs = []
    for strand, direction in ((dna_seq, 1), (reverse_complement(dna_seq), -1)):
        for i in range(len(strand)):
            if strand[i:i+3] in start_codons:
                orf = ""
                start_index = i
                for j in range(i, len(strand), 3):
                    codon = strand[j:j+3]
                    if codon in stop_codons:
                        orf += codon
                        end_index = j + 3
                        if direction == 1:
                            orfs.append((start_index, end_index, direction, orf))
                        else:
                            orfs.append((len(strand) - end_index, len(strand) - start_index, direction, orf))
                        break
                    orf += codon

    # Filter out overlapping ORFs
    filtered_orfs = set()
    for i, orf1 in enumerate(orfs):
        overlap = False
        for j, orf2 in enumerate(orfs):
            if i != j and \
            orf1[0] < orf2[1] and \
            orf1[1] > orf2[0] and \
            orf1[2] == orf2[2]:
                overlap = True
                if len(orf1[3]) > len(orf2[3]):
                    filtered_orfs.add(orf1)
                else:
                    filtered_orfs.add(orf2)
                break
        if not overlap:
            filtered_orfs.add(orf1)

    return sorted(filtered_orfs, key=lambda x: (x[0], x[1], x[2], len(x[3])))

def codon_bias(dna_seq: str) -> dict:
    # Dictionary to store codon frequencies
    codon_counts = {}
    
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1
    
    total_codons = sum(codon_counts.values())
    
    expected_frequency = 1 / (total_codons / 64)  # assuming uniform usage
    
    codon_bias_values = {}
    for codon, count in codon_counts.items():
        expected_count = total_codons * expected_frequency
        codon_bias_values[codon] = count / expected_count
    
    return codon_bias_values

dna_sequence = Dna(id='ID1', sequence='AGCTAGCTAATGTACCATTACATCGTAGCTAGCTAGATATCCTAGCGCGCT')
longest_orfs = find_orf(dna_sequence.sequence)
translated_orfs = [translate(i[3]) for i in longest_orfs]
print(translated_orfs)
