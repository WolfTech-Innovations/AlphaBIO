import curses
import random

# Branding for the program
def print_branding():
    print("AlphaBIO 2.1.3 - by WolfTech Innovations\n")
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

    valid_bases = 'ATCG' if sequence_type == "DNA" else 'AUCGA'

    while not all(base in valid_bases for base in target_sequence):
        stdscr.addstr(f"Invalid input. Please enter a valid {sequence_type} sequence using only {valid_bases}.\n")
        stdscr.refresh()
        target_sequence = stdscr.getstr().decode("utf-8").upper()

    return target_sequence

# Save results to a file
def save_results_to_file(sequence, filename="output.txt"):
    with open(filename, 'w') as file:
        file.write(sequence)

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
        stdscr.addstr("Enter your choice (1/2/3/4/5/6/7/8/9): ")

        stdscr.refresh()
        key = stdscr.getch()

        if key == 49:  # '1' key pressed
            action_choice = 1
        elif key == 50:  # '2' key pressed
            action_choice = 2
        elif key == 51:  # '3' key pressed
            action_choice = 3
        elif key == 52:  # '4' key pressed
            action_choice = 4
        elif key == 53:  # '5' key pressed
            action_choice = 5
        elif key == 54:  # '6' key pressed
            action_choice = 6
        elif key == 55:  # '7' key pressed
            action_choice = 7
        elif key == 56:  # '8' key pressed
            action_choice = 8
        elif key == 57:  # '9' key pressed
            action_choice = 9

    if action_choice == 1 and sequence_type_choice == "DNA":
        rna_sequence = dna_to_rna(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"Transcribed RNA Sequence: {rna_sequence}\n")
    elif action_choice == 2 and sequence_type_choice == "RNA":
        dna_sequence = rna_to_dna(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"Transcribed DNA Sequence: {dna_sequence}\n")
    elif action_choice == 3:
        stdscr.addstr("Choose mutation type: point, insertion, or deletion: ")
        stdscr.refresh()
        mutation_type = stdscr.getstr().decode("utf-8").lower()
        mutated_sequence = introduce_mutation(target_sequence, mutation_type)
        stdscr.clear()
        stdscr.addstr(f"Mutated Sequence: {mutated_sequence}\n")
    elif action_choice == 4:
        reverse_comp = reverse_complement(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"Reverse Complement: {reverse_comp}\n")
    elif action_choice == 5:
        gc_content = calculate_gc_content(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"GC Content: {gc_content:.2f}%\n")
    elif action_choice == 6:
        base_counts = count_bases(target_sequence)
        stdscr.clear()
        stdscr.addstr(f"Base Counts: A={base_counts['A']}, T={base_counts['T']}, C={base_counts['C']}, G={base_counts['G']}\n")
    elif action_choice == 7:
        stdscr.addstr("Enter the desired length for the random sequence: ")
        stdscr.refresh()
        length = int(stdscr.getstr().decode("utf-8"))
        random_sequence = generate_random_sequence(length)
        stdscr.clear()
        stdscr.addstr(f"Generated Random Sequence: {random_sequence}\n")
    elif action_choice == 8:
        stdscr.addstr("Enter a second sequence to align with:\n")
        stdscr.refresh()
        second_sequence = stdscr.getstr().decode("utf-8").upper()
        matches, mismatches = align_sequences(target_sequence, second_sequence)
        stdscr.clear()
        stdscr.addstr(f"Alignment: Matches={matches}, Mismatches={mismatches}\n")
    elif action_choice == 9:
        save_results_to_file(target_sequence, "output.txt")
        stdscr.clear()
        stdscr.addstr("Results saved to 'output.txt'.\n")
    
    stdscr.refresh()
    stdscr.getch()  # Wait for user to press a key before exiting

if __name__ == "__main__":
    curses.wrapper(main)
