import curses
import random

# Branding for the program
def print_branding():
    print("AlphaBIO 2.2.0 - by WolfTech Innovations\n")
    print("*****************************************\n")

# Function to generate complementary DNA strand
def generate_complementary_dna_strand(target_sequence):
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    if not all(base in complement for base in target_sequence):
        raise ValueError("Invalid DNA sequence. Only 'A', 'T', 'C', 'G' are allowed.")
    
    complementary_strand = ''.join([complement[base] for base in target_sequence])
    return complementary_strand

# Function to generate complementary RNA strand
def generate_complementary_rna_strand(target_sequence):
    complement = {
        'A': 'U',
        'U': 'A',
        'C': 'G',
        'G': 'C'
    }

    if not all(base in complement for base in target_sequence):
        raise ValueError("Invalid RNA sequence. Only 'A', 'U', 'C', 'G' are allowed.")

    complementary_strand = ''.join([complement[base] for base in target_sequence])
    return complementary_strand

# Function to transcribe DNA to RNA
def dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')

# Function to transcribe RNA to DNA
def rna_to_dna(rna_sequence):
    return rna_sequence.replace('U', 'T')

# Function to create cDNA from DNA
def dna_to_cdna(dna_sequence):
    # cDNA is the complement of the DNA template strand
    return generate_complementary_dna_strand(dna_sequence)

# Function to create viral treatment sequence
def create_virus_treatment(viral_rna):
    """
    Creates a potential treatment sequence for viral RNA by:
    1. Converting viral RNA to DNA
    2. Creating cDNA from the DNA
    3. Converting cDNA back to RNA (antisense RNA)
    
    This antisense RNA can potentially bind to viral RNA and prevent translation.
    """
    # Step 1: Convert viral RNA to DNA
    dna_sequence = rna_to_dna(viral_rna)
    
    # Step 2: Create cDNA from DNA
    cdna_sequence = dna_to_cdna(dna_sequence)
    
    # Step 3: Convert cDNA back to RNA (antisense RNA)
    antisense_rna = dna_to_rna(cdna_sequence)
    
    return {
        "original_viral_rna": viral_rna,
        "dna_sequence": dna_sequence,
        "cdna_sequence": cdna_sequence,
        "antisense_rna": antisense_rna
    }

# Function to introduce mutations (point mutation, insertion, or deletion)
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

# Function to calculate GC content percentage
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    return (gc_count / total_bases) * 100 if total_bases > 0 else 0

# Function to count each base (A, T/U, C, G)
def count_bases(sequence):
    base_counts = {'A': sequence.count('A'), 'T': sequence.count('T'), 'C': sequence.count('C'), 'G': sequence.count('G')}
    # Add U count for RNA sequences
    if 'U' in sequence:
        base_counts['U'] = sequence.count('U')
        base_counts.pop('T', None)  # Remove T count for RNA
    return base_counts

# Function to reverse complement a sequence
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'} if 'U' not in sequence else {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    reversed_sequence = sequence[::-1]
    return ''.join([complement[base] for base in reversed_sequence])

# Function to get a random sequence of a specific length and GC content
def generate_random_sequence(length, gc_content=50):
    sequence = []
    gc_bases = int(length * gc_content / 100)
    at_bases = length - gc_bases
    
    sequence.extend(['G', 'C'] * (gc_bases // 2))
    sequence.extend(['A', 'T'] * (at_bases // 2))

    random.shuffle(sequence)
    return ''.join(sequence)

# Function to get input from the user for sequences
def get_input_sequence(sequence_type, stdscr):
    if sequence_type == "DNA":
        stdscr.addstr("Enter a DNA sequence (only 'A', 'T', 'C', 'G' allowed):\n")
    elif sequence_type == "RNA":
        stdscr.addstr("Enter an RNA sequence (only 'A', 'U', 'C', 'G' allowed):\n")
    
    stdscr.refresh()
    target_sequence = stdscr.getstr().decode("utf-8").upper()

    valid_bases = 'ATCG' if sequence_type == "DNA" else 'AUCG'

    while not all(base in valid_bases for base in target_sequence):
        stdscr.addstr(f"Invalid input. Please enter a valid {sequence_type} sequence using only {valid_bases}.\n")
        stdscr.refresh()
        target_sequence = stdscr.getstr().decode("utf-8").upper()

    return target_sequence

# Save results to a file
def save_results_to_file(content, filename="output.txt"):
    with open(filename, 'w') as file:
        file.write(content)

# Sequence alignment (simple)
def align_sequences(seq1, seq2):
    # Basic sequence comparison
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return matches, mismatches

# Main Program (with curses-based TUI)
def main(stdscr):
    # Initialize the screen
    curses.curs_set(0)  # Hide the cursor
    stdscr.clear()  # Clear the screen
    print_branding()
    
    # Sequence Type Selector (TUI)
    sequence_type_choice = None
    options = ["DNA", "RNA"]
    
    while sequence_type_choice is None:
        stdscr.clear()
        stdscr.addstr("Select the type of sequence you want to work with:\n")
        for idx, option in enumerate(options):
            stdscr.addstr(f"{idx+1}. {option}\n")
        stdscr.addstr("Enter 1 for DNA or 2 for RNA: ")

        stdscr.refresh()
        key = stdscr.getch()

        if key == 49:  # '1' key pressed
            sequence_type_choice = "DNA"
        elif key == 50:  # '2' key pressed
            sequence_type_choice = "RNA"

    # Get the sequence from the user
    target_sequence = get_input_sequence(sequence_type_choice, stdscr)

    # Generate complementary strand
    if sequence_type_choice == "DNA":
        complementary_strand = generate_complementary_dna_strand(target_sequence)
    else:
        complementary_strand = generate_complementary_rna_strand(target_sequence)

    # TUI Menu for additional actions
    action_choice = None
    actions = [
        "Transcribe DNA to RNA (if DNA input)",
        "Transcribe RNA to DNA (if RNA input)",
        "Create Virus Treatment (RNA to DNA to cDNA to RNA)",
        "Introduce Mutation (point, insertion, or deletion)",
        "Reverse Complement of Sequence",
        "Calculate GC Content",
        "Count Base Pairs (A, T/U, C, G)",
        "Generate Random Sequence",
        "Align Two Sequences",
        "Save Results to File",
        "Exit"
    ]
    
    while action_choice is None:
        stdscr.clear()
        stdscr.addstr("Select an additional action:\n")
        for idx, action in enumerate(actions):
            stdscr.addstr(f"{idx+1}. {action}\n")
        stdscr.addstr("Enter your choice (1-11): ")

        stdscr.refresh()
        key = stdscr.getch()

        if 49 <= key <= 57:  # '1' to '9' keys pressed
            action_choice = key - 48
        elif key == 49 + 10:  # '10' key combination
            action_choice = 10
        elif key == 49 + 11:  # '11' key combination
            action_choice = 11

    if action_choice == 1 and sequence_type_choice == "DNA":
        rna_sequence = dna_to_rna(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"Transcribed RNA Sequence: {rna_sequence}\n")
    elif action_choice == 2 and sequence_type_choice == "RNA":
        dna_sequence = rna_to_dna(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"Transcribed DNA Sequence: {dna_sequence}\n")
    elif action_choice == 3:
        if sequence_type_choice == "RNA":
            treatment = create_virus_treatment(target_sequence)
            stdscr.clear()
            stdscr.addstr("Virus Treatment Created:\n")
            stdscr.addstr(f"Original Viral RNA: {treatment['original_viral_rna']}\n")
            stdscr.addstr(f"DNA Sequence: {treatment['dna_sequence']}\n")
            stdscr.addstr(f"cDNA Sequence: {treatment['cdna_sequence']}\n")
            stdscr.addstr(f"Antisense RNA Treatment: {treatment['antisense_rna']}\n")
            stdscr.addstr("\nThis antisense RNA can bind to viral RNA and potentially inhibit viral replication.\n")
        else:
            stdscr.clear()
            stdscr.addstr("Virus treatments can only be created from RNA sequences.\n")
    elif action_choice == 4:
        stdscr.clear()
        stdscr.addstr("Choose mutation type (point/insertion/deletion): ")
        stdscr.refresh()
        mutation_type = stdscr.getstr().decode("utf-8").lower()
        mutated_sequence = introduce_mutation(target_sequence, mutation_type)
        stdscr.clear()
        stdscr.addstr(f"Mutated Sequence: {mutated_sequence}\n")
    elif action_choice == 5:
        reverse_comp = reverse_complement(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"Reverse Complement: {reverse_comp}\n")
    elif action_choice == 6:
        gc_content = calculate_gc_content(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"GC Content: {gc_content:.2f}%\n")
    elif action_choice == 7:
        base_counts = count_bases(target_sequence)
        stdscr.clear()
        stdscr.addstr("Base Counts:\n")
        for base, count in base_counts.items():
            stdscr.addstr(f"{base}={count} ")
        stdscr.addstr("\n")
    elif action_choice == 8:
        stdscr.clear()
        stdscr.addstr("Enter the desired length for the random sequence: ")
        stdscr.refresh()
        length = int(stdscr.getstr().decode("utf-8"))
        random_sequence = generate_random_sequence(length)
        stdscr.clear()
        stdscr.addstr(f"Generated Random Sequence: {random_sequence}\n")
    elif action_choice == 9:
        stdscr.clear()
        stdscr.addstr("Enter a second sequence to align with:\n")
        stdscr.refresh()
        second_sequence = stdscr.getstr().decode("utf-8").upper()
        matches, mismatches = align_sequences(target_sequence, second_sequence)
        stdscr.clear()
        stdscr.addstr(f"Alignment: Matches={matches}, Mismatches={mismatches}\n")
    elif action_choice == 10:
        if sequence_type_choice == "RNA" and action_choice == 3:
            treatment = create_virus_treatment(target_sequence)
            save_content = f"""Virus Treatment Report
Original Viral RNA: {treatment['original_viral_rna']}
DNA Sequence: {treatment['dna_sequence']}
cDNA Sequence: {treatment['cdna_sequence']}
Antisense RNA Treatment: {treatment['antisense_rna']}

This antisense RNA can bind to viral RNA and potentially inhibit viral replication.
"""
            save_results_to_file(save_content, "viral_treatment.txt")
            stdscr.clear()
            stdscr.addstr("Results saved to 'viral_treatment.txt'.\n")
        else:
            save_results_to_file(target_sequence, "output.txt")
            stdscr.clear()
            stdscr.addstr("Results saved to 'output.txt'.\n")
    
    stdscr.refresh()
    stdscr.getch()  # Wait for user to press a key before exiting

if __name__ == "__main__":
    curses.wrapper(main)
