import curses
import random

def print_branding():
    print("AlphaBIO 2.1.4 - by WolfTech Innovations\n")
    print("*****************************************\n")

def dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')

def rna_to_dna(rna_sequence):
    return rna_sequence.replace('U', 'T')

def dna_to_cdna(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in dna_sequence)

def cdna_to_rna(cdna_sequence):
    return dna_to_rna(cdna_sequence)

def rna_to_cdna_back_to_rna(rna_sequence):
    dna = rna_to_dna(rna_sequence)
    cdna = dna_to_cdna(dna)
    return cdna_to_rna(cdna)

def generate_complementary_dna_strand(target_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in target_sequence)

def generate_complementary_rna_strand(target_sequence):
    complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in target_sequence)

def introduce_mutation(sequence, mutation_type="point", mutation_probability=0.1):
    mutated_sequence = list(sequence)
    if mutation_type == "point" and random.random() < mutation_probability:
        idx = random.randint(0, len(sequence) - 1)
        new_base = random.choice(['A', 'T', 'C', 'G'])
        mutated_sequence[idx] = new_base
    elif mutation_type == "insertion" and random.random() < mutation_probability:
        idx = random.randint(0, len(sequence))
        new_base = random.choice(['A', 'T', 'C', 'G'])
        mutated_sequence.insert(idx, new_base)
    elif mutation_type == "deletion" and random.random() < mutation_probability:
        idx = random.randint(0, len(sequence) - 1)
        mutated_sequence.pop(idx)
    return ''.join(mutated_sequence)

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    return (gc_count / total_bases) * 100 if total_bases > 0 else 0

def count_bases(sequence):
    return {base: sequence.count(base) for base in "ATCGU" if base in sequence}

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'}
    return ''.join(complement[base] for base in reversed(sequence))

def generate_random_sequence(length, gc_content=50):
    sequence = random.choices(['A', 'T', 'C', 'G'], k=length)
    return ''.join(sequence)

def align_sequences(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return matches, mismatches

def main():
    print_branding()
    sequence_type = input("Enter sequence type (DNA/RNA): ").upper()
    if sequence_type not in ["DNA", "RNA"]:
        print("Invalid sequence type.")
        return
    sequence = input(f"Enter {sequence_type} sequence: ").upper()
    
    if sequence_type == "DNA":
        comp_strand = generate_complementary_dna_strand(sequence)
        print(f"Complementary DNA strand: {comp_strand}")
        print(f"Transcribed RNA: {dna_to_rna(sequence)}")
    else:
        comp_strand = generate_complementary_rna_strand(sequence)
        print(f"Complementary RNA strand: {comp_strand}")
        print(f"Converted DNA: {rna_to_dna(sequence)}")
        print(f"Back to RNA (cDNA to RNA): {rna_to_cdna_back_to_rna(sequence)}")
    
    print(f"Reverse Complement: {reverse_complement(sequence)}")
    print(f"GC Content: {calculate_gc_content(sequence):.2f}%")
    print(f"Base Counts: {count_bases(sequence)}")
    
    mutation = introduce_mutation(sequence)
    print(f"Mutated Sequence: {mutation}")
    
    rand_seq = generate_random_sequence(20)
    print(f"Random Sequence: {rand_seq}")
    
    seq2 = input("Enter another sequence to align: ").upper()
    matches, mismatches = align_sequences(sequence, seq2)
    print(f"Alignment - Matches: {matches}, Mismatches: {mismatches}")

if __name__ == "__main__":
    main()
