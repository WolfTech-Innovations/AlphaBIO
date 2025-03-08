# Branding for the program
def print_branding():
    print("----------------------------------------\n") 
    print("AlphaBIO 2.1.3 - by WolfTech Innovations\n")
    print("----------------------------------------\n")

# Function to generate complementary DNA strand
def generate_complementary_strand(target_sequence):
    # Mapping of bases to their complements
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    # Ensure input sequence is valid
    if not all(base in complement for base in target_sequence):
        raise ValueError("Invalid DNA sequence. Only 'A', 'T', 'C', 'G' are allowed.")

    # Create the complementary strand by mapping each base
    complementary_strand = ''.join([complement[base] for base in target_sequence])
    return complementary_strand

# Function to "disable" (introduce mutations to disrupt a target sequence)
def disable_function(target_sequence):
    import random

    # Randomly introduce a mutation by changing one base at a random position
    disable_index = random.randint(0, len(target_sequence) - 1)
    base_options = ['A', 'T', 'C', 'G']
    
    # Change the base at the randomly selected index to a different one
    new_base = random.choice([base for base in base_options if base != target_sequence[disable_index]])
    disabled_sequence = target_sequence[:disable_index] + new_base + target_sequence[disable_index+1:]
    
    return disabled_sequence

# Function to get realistic data or user input
def get_input_sequence():
    print("Enter a DNA sequence (only 'A', 'T', 'C', 'G' allowed):")
    target_sequence = input().upper()

    # Validate sequence input
    while not all(base in 'ATCG' for base in target_sequence):
        print("Invalid input. Please enter a valid DNA sequence using only 'A', 'T', 'C', 'G'.")
        target_sequence = input().upper()

    return target_sequence

# Main Program
def main():
    print_branding()
    
    # Step 1: Get the target sequence from user input
    target_sequence = get_input_sequence()
    
    # Step 2: Generate the complementary DNA strand
    complementary_strand = generate_complementary_strand(target_sequence)
    
    # Step 3: Introduce a mutation to simulate disabling functionality
    disabled_sequence = disable_function(complementary_strand)
    
    # Output the results
    print("\n********* Results *********")
    print("Target Sequence:       ", target_sequence)
    print("Complementary Strand:  ", complementary_strand)
    print("Disabled (mutated) Sequence: ", disabled_sequence)
    print("\n****************************")
    print("AlphaBIO 2.1.3 - by WolfTech Innovations\n")
    
if __name__ == "__main__":
    main()
