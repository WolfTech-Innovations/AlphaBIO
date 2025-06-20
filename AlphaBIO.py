import curses
import random
import re
import os

# Branding
def print_branding(stdscr):
    stdscr.addstr("AlphaBIO 2.3.0 - by WolfTech Innovations\n")
    stdscr.addstr("**************************************************\n")

# Sequence utilities
def validate_sequence(seq, allowed):
    seq = seq.upper().replace(" ", "")
    if not all(base in allowed for base in seq):
        raise ValueError(f"Invalid sequence: only {allowed} allowed.")
    return seq

def reverse_complement(seq):
    comp = {'A':'T','T':'A','C':'G','G':'C','U':'A'}
    return ''.join(comp[base] for base in seq[::-1])

def generate_complementary_dna_strand(seq):
    return ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[base] for base in seq)

def generate_complementary_rna_strand(seq):
    return ''.join({'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}[base] for base in seq)

def dna_to_rna(seq):
    return seq.replace('T', 'U')

def rna_to_dna(seq):
    return seq.replace('U', 'T')

def dna_to_cdna(seq):
    return generate_complementary_dna_strand(seq)

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    return (gc_count / total_bases) * 100 if total_bases > 0 else 0

def count_bases(sequence):
    base_counts = {'A': sequence.count('A'), 'T': sequence.count('T'), 'C': sequence.count('C'), 'G': sequence.count('G')}
    if 'U' in sequence:
        base_counts['U'] = sequence.count('U')
        base_counts.pop('T', None)
    return base_counts

def introduce_mutation(sequence, mutation_type="point", mutation_probability=0.1):
    mutated_sequence = list(sequence)
    bases = ['A', 'T', 'C', 'G'] if 'U' not in sequence else ['A', 'U', 'C', 'G']
    if mutation_type == "point" and random.random() < mutation_probability:
        idx = random.randint(0, len(sequence) - 1)
        new_base = random.choice([b for b in bases if b != sequence[idx]])
        mutated_sequence[idx] = new_base
    elif mutation_type == "insertion" and random.random() < mutation_probability:
        idx = random.randint(0, len(sequence))
        new_base = random.choice(bases)
        mutated_sequence.insert(idx, new_base)
    elif mutation_type == "deletion" and random.random() < mutation_probability and len(mutated_sequence) > 1:
        idx = random.randint(0, len(sequence) - 1)
        mutated_sequence.pop(idx)
    return ''.join(mutated_sequence)

def generate_random_sequence(length, gc_content=50):
    gc_bases = int(length * gc_content / 100)
    at_bases = length - gc_bases
    sequence = []
    for _ in range(gc_bases):
        sequence.append(random.choice(['G', 'C']))
    for _ in range(at_bases):
        sequence.append(random.choice(['A', 'T']))
    random.shuffle(sequence)
    return ''.join(sequence)

def align_sequences(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    extra = abs(len(seq1) - len(seq2))
    mismatches += extra
    return matches, mismatches

def save_results_to_file(content, filename="output.txt"):
    with open(filename, 'w') as file:
        file.write(content)

# CRISPR functions
def find_pam_sites(dna_seq, pam="NGG"):
    """Returns list of (start, gRNA) tuples for all NGG sites (SpCas9)."""
    dna_seq = dna_seq.upper()
    results = []
    for m in re.finditer(r'(?=([ACGT]{20}GG))', dna_seq):
        gRNA = m.group(1)[:20]
        start = m.start(1)
        results.append((start, gRNA))
    return results

def design_gRNA(dna_seq, target_seq):
    """Return gRNA for target_seq if it has adjacent PAM (NGG)."""
    dna_seq = dna_seq.upper()
    target_seq = target_seq.upper()
    idx = dna_seq.find(target_seq)
    if idx == -1:
        raise ValueError("Target sequence not found in DNA")
    pam_site = idx + len(target_seq)
    if pam_site + 2 <= len(dna_seq) and dna_seq[pam_site:pam_site+2] == "GG":
        gRNA = target_seq[-20:] if len(target_seq) >= 20 else target_seq
        return gRNA, pam_site
    else:
        raise ValueError("No adjacent PAM site found for CRISPR/Cas9 targeting.")

def perform_genome_edit(dna_seq, target_seq, edit_seq):
    """Edit target_seq in dna_seq to edit_seq, and return edited sequence and gRNA."""
    dna_seq = dna_seq.upper()
    target_seq = target_seq.upper()
    edit_seq = edit_seq.upper()
    idx = dna_seq.find(target_seq)
    if idx == -1:
        raise ValueError("Target not found.")
    edited = dna_seq[:idx] + edit_seq + dna_seq[idx + len(target_seq):]
    gRNA, pam_site = design_gRNA(dna_seq, target_seq)
    return edited, gRNA

def create_antisense_crna(rna_seq):
    """Generate antisense cRNA for a given RNA sequence."""
    comp = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp[base] for base in rna_seq[::-1])

def crispr_deliver_antisense(dna_seq, insertion_site, antisense_seq):
    """Simulate CRISPR delivery of antisense cRNA at a cut site."""
    return dna_seq[:insertion_site] + antisense_seq + dna_seq[insertion_site:]

# TUI functions
def get_input_sequence(sequence_type, stdscr):
    allowed = 'ATCG' if sequence_type == "DNA" else 'AUCG'
    prompt = f"Enter a {sequence_type} sequence (only {allowed} allowed):\n"
    while True:
        stdscr.addstr(prompt)
        stdscr.refresh()
        s = stdscr.getstr().decode("utf-8").upper().replace(" ", "")
        try:
            return validate_sequence(s, allowed)
        except ValueError as e:
            stdscr.addstr(str(e) + "\n")
            stdscr.refresh()

def get_integer_input(stdscr, prompt, minval=None, maxval=None):
    while True:
        stdscr.addstr(prompt)
        stdscr.refresh()
        try:
            val = int(stdscr.getstr().decode("utf-8"))
            if (minval is not None and val < minval) or (maxval is not None and val > maxval):
                stdscr.addstr(f"Value must be between {minval} and {maxval}.\n")
            else:
                return val
        except ValueError:
            stdscr.addstr("Please enter a valid integer.\n")

def main(stdscr):
    curses.curs_set(0)
    stdscr.clear()
    print_branding(stdscr)
    stdscr.refresh()

    # Sequence Type Selector
    sequence_type_choice = None
    options = ["DNA", "RNA"]
    while sequence_type_choice is None:
        stdscr.addstr("Select the type of sequence you want to work with:\n")
        for idx, option in enumerate(options):
            stdscr.addstr(f"{idx+1}. {option}\n")
        stdscr.addstr("Enter 1 for DNA or 2 for RNA: ")
        stdscr.refresh()
        key = stdscr.getch()
        if key == 49:  # '1'
            sequence_type_choice = "DNA"
        elif key == 50:  # '2'
            sequence_type_choice = "RNA"
        stdscr.clear()

    target_sequence = get_input_sequence(sequence_type_choice, stdscr)

    # Menu
    actions = [
        "Transcribe DNA to RNA (if DNA input)",
        "Transcribe RNA to DNA (if RNA input)",
        "Create Virus Treatment (RNA to DNA to cDNA to RNA)",
        "Introduce Mutation (point/insertion/deletion)",
        "Reverse Complement of Sequence",
        "Calculate GC Content",
        "Count Base Pairs (A, T/U, C, G)",
        "Generate Random Sequence",
        "Align Two Sequences",
        "Save Results to File",
        "CRISPR: Find PAM sites",
        "CRISPR: Design gRNA for a target",
        "CRISPR: Genome Editing (find/replace sequence)",
        "CRISPR: Deliver cRNA Antisense treatment",
        "Exit"
    ]
    while True:
        stdscr.clear()
        print_branding(stdscr)
        stdscr.addstr("Select an additional action:\n")
        for idx, action in enumerate(actions):
            stdscr.addstr(f"{idx+1}. {action}\n")
        stdscr.addstr(f"Enter your choice (1-{len(actions)}): ")
        stdscr.refresh()
        try:
            action_choice = int(stdscr.getstr().decode("utf-8"))
        except ValueError:
            continue

        stdscr.clear()

        # Standard actions
        if action_choice == 1 and sequence_type_choice == "DNA":
            rna_sequence = dna_to_rna(target_sequence)
            stdscr.addstr(f"Transcribed RNA Sequence: {rna_sequence}\n")
        elif action_choice == 2 and sequence_type_choice == "RNA":
            dna_sequence = rna_to_dna(target_sequence)
            stdscr.addstr(f"Transcribed DNA Sequence: {dna_sequence}\n")
        elif action_choice == 3:
            if sequence_type_choice == "RNA":
                treatment = {}
                treatment['original_viral_rna'] = target_sequence
                treatment['dna_sequence'] = rna_to_dna(target_sequence)
                treatment['cdna_sequence'] = dna_to_cdna(treatment['dna_sequence'])
                treatment['antisense_rna'] = dna_to_rna(treatment['cdna_sequence'])
                stdscr.addstr("Virus Treatment Created:\n")
                for k, v in treatment.items():
                    stdscr.addstr(f"{k}: {v}\n")
                stdscr.addstr("\nThis antisense RNA can bind to viral RNA and potentially inhibit viral replication.\n")
            else:
                stdscr.addstr("Virus treatments can only be created from RNA sequences.\n")
        elif action_choice == 4:
            stdscr.addstr("Choose mutation type (point/insertion/deletion): ")
            stdscr.refresh()
            mutation_type = stdscr.getstr().decode("utf-8").lower()
            mutated_sequence = introduce_mutation(target_sequence, mutation_type)
            stdscr.addstr(f"Mutated Sequence: {mutated_sequence}\n")
        elif action_choice == 5:
            reverse_comp = reverse_complement(target_sequence)
            stdscr.addstr(f"Reverse Complement: {reverse_comp}\n")
        elif action_choice == 6:
            gc_content = calculate_gc_content(target_sequence)
            stdscr.addstr(f"GC Content: {gc_content:.2f}%\n")
        elif action_choice == 7:
            base_counts = count_bases(target_sequence)
            stdscr.addstr("Base Counts:\n")
            for base, count in base_counts.items():
                stdscr.addstr(f"{base}={count} ")
            stdscr.addstr("\n")
        elif action_choice == 8:
            length = get_integer_input(stdscr, "Enter the desired length for the random sequence: ", 1)
            random_sequence = generate_random_sequence(length)
            stdscr.addstr(f"Generated Random Sequence: {random_sequence}\n")
        elif action_choice == 9:
            stdscr.addstr("Enter a second sequence to align with:\n")
            stdscr.refresh()
            second_sequence = stdscr.getstr().decode("utf-8").upper().replace(" ", "")
            matches, mismatches = align_sequences(target_sequence, second_sequence)
            stdscr.addstr(f"Alignment: Matches={matches}, Mismatches={mismatches}\n")
        elif action_choice == 10:
            filename = "output.txt"
            save_results_to_file(target_sequence, filename)
            stdscr.addstr(f"Results saved to '{filename}'.\n")
        # CRISPR: Find PAM sites
        elif action_choice == 11 and sequence_type_choice == "DNA":
            pam_sites = find_pam_sites(target_sequence)
            if pam_sites:
                stdscr.addstr("PAM sites (NGG) and their gRNAs:\n")
                for idx, (start, gRNA) in enumerate(pam_sites):
                    stdscr.addstr(f"{idx+1}: Start={start}, gRNA={gRNA}\n")
            else:
                stdscr.addstr("No PAM (NGG) sites found in the DNA.\n")
        # CRISPR: Design gRNA for a target
        elif action_choice == 12 and sequence_type_choice == "DNA":
            stdscr.addstr("Enter a target subsequence (must be in DNA):\n")
            stdscr.refresh()
            target_subseq = stdscr.getstr().decode("utf-8").upper().replace(" ", "")
            try:
                gRNA, pam_site = design_gRNA(target_sequence, target_subseq)
                stdscr.addstr(f"gRNA: {gRNA}, PAM at position {pam_site}\n")
            except Exception as e:
                stdscr.addstr(f"Error: {e}\n")
        # CRISPR: Genome Editing (find/replace sequence)
        elif action_choice == 13 and sequence_type_choice == "DNA":
            stdscr.addstr("Enter the subsequence to REPLACE (must be present):\n")
            stdscr.refresh()
            old_seq = stdscr.getstr().decode("utf-8").upper().replace(" ", "")
            stdscr.addstr("Enter the NEW sequence:\n")
            stdscr.refresh()
            new_seq = stdscr.getstr().decode("utf-8").upper().replace(" ", "")
            try:
                edited, gRNA = perform_genome_edit(target_sequence, old_seq, new_seq)
                stdscr.addstr(f"Edited DNA: {edited}\ngRNA: {gRNA}\n")
            except Exception as e:
                stdscr.addstr(f"Error: {e}\n")
        # CRISPR: Deliver cRNA Antisense treatment
        elif action_choice == 14:
            if sequence_type_choice == "RNA":
                antisense = create_antisense_crna(target_sequence)
                stdscr.addstr(f"Antisense cRNA: {antisense}\n")
                if sequence_type_choice == "DNA":
                    stdscr.addstr("Enter insertion site for antisense (0-based index): ")
                    stdscr.refresh()
                    site = get_integer_input(stdscr, "", 0, len(target_sequence))
                    edited_genome = crispr_deliver_antisense(target_sequence, site, antisense)
                    stdscr.addstr(f"Genome after cRNA Antisense Treatment: {edited_genome}\n")
                else:
                    stdscr.addstr("Antisense RNA generated (CRISPR delivery simulated only for DNA).\n")
            else:
                stdscr.addstr("Antisense cRNA can only be generated from RNA sequences.\n")
        elif action_choice == 15:
            stdscr.addstr("Exiting...\n")
            stdscr.refresh()
            curses.napms(1000)
            break
        else:
            stdscr.addstr("Invalid choice or action not available for this sequence type.\n")
        stdscr.addstr("\nPress any key to continue...")
        stdscr.refresh()
        stdscr.getch()

if __name__ == "__main__":
    curses.wrapper(main)
