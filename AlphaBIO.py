"""
AlphaBIO™ DNA Binding Designer
Version 1.2.3
© 2025 WolfTech Innovations
"""

import re
import numpy as np
from datetime import datetime

class AlphaBIODesigner:
    """
    AlphaBIO™ DNA Binding Designer by WolfTech Innovations
    
    A professional tool for designing nucleic acid-based therapeutics that can bind to and 
    inactivate target DNA sequences using realistic biological strategies.
    """
    
    def __init__(self):
        # Initialize with valid nucleotides and complementary pairs
        self.valid_bases = {'A', 'T', 'G', 'C'}
        self.complementary_pairs = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G'
        }
        
        # Binding parameters based on real biochemical properties
        self.gc_content_optimal = (0.45, 0.65)  # Optimal GC content range
        self.melting_temp_optimal = (60, 75)    # Optimal melting temperature range (°C)
        self.self_complementarity_threshold = 6  # Maximum self-complementary bases allowed
        
        # Print startup banner
        self._print_banner()
    
    def _print_banner(self):
        """Display the AlphaBIO startup banner."""
        banner = """
        ╔═════════════════════════════════════════════════╗
        ║                                                 ║
        ║         AlphaBIO™ DNA Binding Designer         ║
        ║             WolfTech Innovations               ║
        ║                 Version 1.2.3                  ║
        ║                                                 ║
        ╚═════════════════════════════════════════════════╝
        """
        print(banner)
        print(f"Session started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("Ready for sequence analysis and binding design.\n")
    
    def validate_sequence(self, sequence):
        """Validate that the input sequence contains only valid DNA bases."""
        if not sequence:
            return False, "Sequence cannot be empty."
        
        sequence = sequence.upper().replace(" ", "")
        
        # Check for valid DNA characters
        if not all(base in self.valid_bases for base in sequence):
            invalid_chars = set(sequence) - self.valid_bases
            return False, f"Invalid nucleotides detected: {', '.join(invalid_chars)}"
        
        return True, sequence
    
    def calculate_gc_content(self, sequence):
        """Calculate the GC content of a sequence."""
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence) if len(sequence) > 0 else 0
    
    def calculate_melting_temp(self, sequence):
        """
        Calculate the approximate melting temperature using the nearest-neighbor method.
        This is a simplified version of the actual thermodynamic calculations.
        """
        sequence = sequence.upper()
        length = len(sequence)
        
        if length < 14:
            # Wallace rule for short oligos
            return 2 * (sequence.count('A') + sequence.count('T')) + 4 * (sequence.count('G') + sequence.count('C'))
        else:
            # Simplified nearest-neighbor method for longer sequences
            gc_content = self.calculate_gc_content(sequence)
            return 81.5 + 16.6 * np.log10(0.05) + 41 * gc_content - (600 / length)
    
    def check_self_complementarity(self, sequence):
        """Check for self-complementary regions within a sequence."""
        sequence = sequence.upper()
        comp_sequence = self.generate_complementary_strand(sequence)
        
        max_complementary = 0
        for i in range(len(sequence) - 3):  # Minimum 4 bases for meaningful complementarity
            for j in range(4, min(len(sequence) - i + 1, 10)):  # Check segments of 4-10 bases
                segment = sequence[i:i+j]
                rev_comp = ''.join([self.complementary_pairs[b] for b in segment[::-1]])
                
                if rev_comp in sequence:
                    max_complementary = max(max_complementary, j)
        
        return max_complementary
    
    def generate_complementary_strand(self, sequence):
        """Generate the complementary DNA strand for the given sequence."""
        return ''.join(self.complementary_pairs[base] for base in sequence.upper())
    
    def design_antisense_oligos(self, target_sequence, min_length=18, max_length=25):
        """
        Design antisense oligonucleotides (ASOs) that can bind to the target mRNA sequence.
        ASOs are short synthetic single-stranded DNA that can bind to mRNA and prevent translation.
        """
        valid, sequence = self.validate_sequence(target_sequence)
        if not valid:
            return {"error": sequence}, "Invalid sequence"
        
        # Find optimal regions for targeting based on:
        # 1. Start codon region (translation initiation)
        # 2. Regions with optimal predicted accessibility
        # 3. Avoid regions with strong secondary structures
        
        candidates = []
        sequence_length = len(sequence)
        
        # Scan the sequence with sliding windows of different lengths
        for length in range(min_length, min(max_length + 1, sequence_length + 1)):
            for i in range(0, sequence_length - length + 1):
                potential_target = sequence[i:i+length]
                
                # Generate the antisense oligonucleotide
                antisense_oligo = self.generate_complementary_strand(potential_target)
                
                # Calculate properties
                gc_content = self.calculate_gc_content(antisense_oligo)
                melting_temp = self.calculate_melting_temp(antisense_oligo)
                self_comp = self.check_self_complementarity(antisense_oligo)
                
                # Score the candidate based on its properties
                score = self._score_oligo(
                    antisense_oligo, 
                    gc_content, 
                    melting_temp, 
                    self_comp,
                    i,  # Position in the target sequence
                    sequence_length
                )
                
                candidates.append({
                    "target_region": (i, i+length-1),
                    "target_sequence": potential_target,
                    "antisense_oligo": antisense_oligo,
                    "length": length,
                    "gc_content": gc_content,
                    "melting_temp": melting_temp,
                    "self_complementarity": self_comp,
                    "score": score
                })
        
        # Sort candidates by score (higher is better)
        candidates.sort(key=lambda x: x["score"], reverse=True)
        
        # Take the top 5 candidates
        top_candidates = candidates[:5]
        
        return {
            "strategy": "antisense oligonucleotides",
            "target_sequence": sequence,
            "target_length": sequence_length,
            "candidates": top_candidates,
            "design_notes": "Antisense oligonucleotides designed to bind to target sequence and inhibit translation"
        }, "Successfully designed antisense oligonucleotides"
    
    def _score_oligo(self, oligo, gc_content, melting_temp, self_comp, position, seq_length):
        """Score an oligonucleotide based on its properties and position."""
        score = 100  # Start with a perfect score
        
        # Penalty for suboptimal GC content
        if gc_content < self.gc_content_optimal[0]:
            score -= 10 * (self.gc_content_optimal[0] - gc_content)
        elif gc_content > self.gc_content_optimal[1]:
            score -= 10 * (gc_content - self.gc_content_optimal[1])
        
        # Penalty for suboptimal melting temperature
        if melting_temp < self.melting_temp_optimal[0]:
            score -= 2 * (self.melting_temp_optimal[0] - melting_temp)
        elif melting_temp > self.melting_temp_optimal[1]:
            score -= 2 * (melting_temp - self.melting_temp_optimal[1])
        
        # Penalty for self-complementarity
        if self_comp > self.self_complementarity_threshold:
            score -= 5 * (self_comp - self.self_complementarity_threshold)
        
        # Bonus for targeting 5' end (start codon region)
        if position < seq_length / 3:
            score += 15 * (1 - position / (seq_length / 3))
        
        # Penalty for repeating bases (e.g., AAAA or GGGG)
        for base in self.valid_bases:
            repeat = base * 4  # Four or more repeats
            if repeat in oligo:
                score -= 10
        
        return max(0, score)  # Score cannot be negative
    
    def design_triplex_forming_oligos(self, target_sequence, min_length=15, max_length=25):
        """
        Design triplex-forming oligonucleotides (TFOs) that can bind to the major groove of DNA.
        TFOs form Hoogsteen or reverse Hoogsteen hydrogen bonds with purine-rich regions of dsDNA.
        """
        valid, sequence = self.validate_sequence(target_sequence)
        if not valid:
            return {"error": sequence}, "Invalid sequence"
        
        # TFOs typically target homopurine/homopyrimidine stretches
        # Find purine-rich regions (A and G)
        purine_rich_regions = []
        
        for i in range(len(sequence) - min_length + 1):
            window = sequence[i:i+min_length]
            purine_count = window.count('A') + window.count('G')
            
            if purine_count >= 0.65 * min_length:  # At least 65% purines
                for j in range(min_length, min(max_length + 1, len(sequence) - i + 1)):
                    extended_window = sequence[i:i+j]
                    extended_purine_count = extended_window.count('A') + extended_window.count('G')
                    
                    if extended_purine_count >= 0.65 * j:
                        purine_rich_regions.append((i, i+j-1, extended_window))
        
        if not purine_rich_regions:
            return {
                "strategy": "triplex-forming oligonucleotides",
                "target_sequence": sequence,
                "candidates": [],
                "design_notes": "No suitable purine-rich regions found for TFO design"
            }, "No suitable regions for TFO design"
        
        # Design TFOs for each suitable region
        tfo_candidates = []
        
        for start, end, region in purine_rich_regions:
            # Design a parallel TFO (binds in parallel orientation)
            # Rules: C binds to G, T binds to A
            parallel_tfo = ''
            for base in region:
                if base == 'G':
                    parallel_tfo += 'C'
                elif base == 'A':
                    parallel_tfo += 'T'
                else:
                    parallel_tfo += 'X'  # X represents a position that doesn't form a triplex
            
            # Calculate properties
            gc_content = self.calculate_gc_content(parallel_tfo.replace('X', ''))
            melting_temp = self.calculate_melting_temp(parallel_tfo.replace('X', ''))
            
            # Score based on length of usable sequence (without X)
            usable_length = len(parallel_tfo) - parallel_tfo.count('X')
            score = 100 * (usable_length / len(parallel_tfo))
            
            tfo_candidates.append({
                "target_region": (start, end),
                "target_sequence": region,
                "tfo_sequence": parallel_tfo,
                "usable_length": usable_length,
                "total_length": len(parallel_tfo),
                "gc_content": gc_content,
                "melting_temp": melting_temp,
                "score": score
            })
        
        # Sort candidates by score
        tfo_candidates.sort(key=lambda x: x["score"], reverse=True)
        
        return {
            "strategy": "triplex-forming oligonucleotides",
            "target_sequence": sequence,
            "target_length": len(sequence),
            "candidates": tfo_candidates[:5],  # Top 5 candidates
            "design_notes": "TFOs designed to form triplex structures with purine-rich regions of the target DNA"
        }, "Successfully designed triplex-forming oligonucleotides"
    
    def design_crispr_guide_rnas(self, target_sequence, pam_sequence="NGG"):
        """
        Design CRISPR guide RNAs to target the specified DNA sequence.
        Default PAM sequence is NGG for SpCas9.
        """
        valid, sequence = self.validate_sequence(target_sequence)
        if not valid:
            return {"error": sequence}, "Invalid sequence"
        
        # Define parameters
        guide_length = 20  # Standard gRNA length for SpCas9
        
        # Find all potential PAM sites
        pam_sites = []
        pam_pattern = pam_sequence.replace('N', '.')
        
        for match in re.finditer(pam_pattern, sequence):
            pam_pos = match.start()
            
            # Guide RNA targets the 20 nucleotides upstream of the PAM
            if pam_pos >= guide_length:
                target_start = pam_pos - guide_length
                target_end = pam_pos - 1
                target_seq = sequence[target_start:pam_pos]
                
                # Calculate properties
                gc_content = self.calculate_gc_content(target_seq)
                
                # Check for poly-T (termination signal for RNA Pol III)
                has_poly_t = 'TTTT' in target_seq
                
                # Simple off-target risk assessment (avoid guides with extreme GC content)
                off_target_risk = "High" if gc_content < 0.3 or gc_content > 0.7 else "Medium" if gc_content < 0.4 or gc_content > 0.6 else "Low"
                
                # Score the guide RNA
                score = self._score_guide_rna(target_seq, gc_content, has_poly_t)
                
                pam_sites.append({
                    "target_region": (target_start, target_end),
                    "target_sequence": target_seq,
                    "pam_position": pam_pos,
                    "pam_sequence": sequence[pam_pos:pam_pos+len(pam_sequence)],
                    "full_target": f"{target_seq}+{sequence[pam_pos:pam_pos+len(pam_sequence)]}",
                    "gc_content": gc_content,
                    "has_poly_t": has_poly_t,
                    "off_target_risk": off_target_risk,
                    "score": score
                })
        
        if not pam_sites:
            return {
                "strategy": "CRISPR guide RNAs",
                "target_sequence": sequence,
                "candidates": [],
                "design_notes": f"No suitable PAM sites ({pam_sequence}) found in the target sequence"
            }, "No suitable PAM sites found"
        
        # Sort by score
        pam_sites.sort(key=lambda x: x["score"], reverse=True)
        
        return {
            "strategy": "CRISPR guide RNAs",
            "target_sequence": sequence,
            "target_length": len(sequence),
            "pam_sequence": pam_sequence,
            "candidates": pam_sites[:5],  # Top 5 candidates
            "design_notes": "Guide RNAs designed to direct Cas9 to specific regions of the target DNA"
        }, "Successfully designed CRISPR guide RNAs"
    
    def _score_guide_rna(self, guide_seq, gc_content, has_poly_t):
        """Score a guide RNA based on its properties."""
        score = 100
        
        # Penalty for extreme GC content
        if gc_content < 0.3:
            score -= 30 * (0.3 - gc_content)
        elif gc_content > 0.7:
            score -= 30 * (gc_content - 0.7)
        
        # Major penalty for poly-T (causes premature termination)
        if has_poly_t:
            score -= 50
        
        # Penalty for G at position 1 (can affect U6 transcription)
        if guide_seq[0] == 'G':
            score -= 5
        
        # Preference for G/C in the last 4 positions (seed region stability)
        seed_region = guide_seq[-4:]
        seed_gc = seed_region.count('G') + seed_region.count('C')
        if seed_gc < 2:
            score -= 10
        
        return max(0, score)
    
    def generate_report(self, target_sequence):
        """Generate a comprehensive binding strategy report for the target sequence."""
        valid, sequence = self.validate_sequence(target_sequence)
        if not valid:
            return {"error": sequence}, "Invalid sequence"
        
        # Analyze sequence properties
        sequence_length = len(sequence)
        gc_content = self.calculate_gc_content(sequence)
        
        # Design different binding strategies
        antisense_results, _ = self.design_antisense_oligos(sequence)
        triplex_results, _ = self.design_triplex_forming_oligos(sequence)
        crispr_results, _ = self.design_crispr_guide_rnas(sequence)
        
        # Determine best overall strategy based on candidates
        strategies = [
            (antisense_results, "antisense oligonucleotides", 
             len(antisense_results.get("candidates", [])) and antisense_results.get("candidates", [{}])[0].get("score", 0)),
            (triplex_results, "triplex-forming oligonucleotides", 
             len(triplex_results.get("candidates", [])) and triplex_results.get("candidates", [{}])[0].get("score", 0)),
            (crispr_results, "CRISPR guide RNAs", 
             len(crispr_results.get("candidates", [])) and crispr_results.get("candidates", [{}])[0].get("score", 0))
        ]
        
        # Sort strategies by score
        strategies.sort(key=lambda x: x[2], reverse=True)
        
        # Compile the report
        report = {
            "analysis_timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "software": "AlphaBIO™ DNA Binding Designer",
            "version": "1.2.3",
            "company": "WolfTech Innovations",
            
            "target_sequence": sequence,
            "sequence_length": sequence_length,
            "gc_content": gc_content,
            
            "recommended_strategy": strategies[0][1] if strategies[0][2] > 0 else "No viable strategy found",
            "strategy_score": strategies[0][2] if strategies[0][2] > 0 else 0,
            
            "all_strategies": {
                "antisense_oligonucleotides": antisense_results,
                "triplex_forming_oligonucleotides": triplex_results,
                "crispr_guide_rnas": crispr_results
            },
            
            "conclusion": f"Based on sequence analysis, the recommended approach for targeting this DNA sequence is using {strategies[0][1]}."
        }
        
        return report, f"Generated comprehensive binding strategy report for {sequence_length}bp sequence"


# Example usage
def main():
    designer = AlphaBIODesigner()
    
    print("Welcome to AlphaBIO™ DNA Binding Designer by WolfTech Innovations")
    print("This tool designs nucleic acid-based therapeutics to target specific DNA sequences.")
    
    while True:
        print("\nPlease enter a DNA sequence to analyze (or 'exit' to quit):")
        user_input = input("> ")
        
        if user_input.lower() == 'exit':
            print("\nThank you for using AlphaBIO™ DNA Binding Designer.")
            print("© 2025 WolfTech Innovations. All rights reserved.")
            break
        
        valid, message = designer.validate_sequence(user_input)
        if valid:
            print(f"\nAnalyzing sequence ({len(message)} bp)...")
            report, status = designer.generate_report(message)
            
            print("\n" + "="*80)
            print(f"ANALYSIS REPORT - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print("="*80)
            
            print(f"Target sequence: {report['target_sequence'][:20]}...{report['target_sequence'][-20:] if len(report['target_sequence']) > 40 else ''}")
            print(f"Sequence length: {report['sequence_length']} bp")
            print(f"GC content: {report['gc_content']:.2f}")
            print(f"\nRECOMMENDED STRATEGY: {report['recommended_strategy']} (Score: {report['strategy_score']:.1f})")
            
            # Display top candidate from recommended strategy
            recommended = report['recommended_strategy'].replace(" ", "_")
            if recommended in ["antisense_oligonucleotides", "triplex_forming_oligonucleotides", "crispr_guide_rnas"]:
                strategy_data = report['all_strategies'][recommended]
                if strategy_data.get("candidates") and len(strategy_data["candidates"]) > 0:
                    top_candidate = strategy_data["candidates"][0]
                    print("\nTop binding candidate:")
                    
                    if "antisense_oligo" in top_candidate:
                        print(f"Antisense oligo: {top_candidate['antisense_oligo']}")
                        print(f"Target region: {top_candidate['target_region']}")
                        print(f"GC content: {top_candidate['gc_content']:.2f}")
                        print(f"Melting temp: {top_candidate['melting_temp']:.1f}°C")
                    
                    elif "tfo_sequence" in top_candidate:
                        print(f"TFO sequence: {top_candidate['tfo_sequence']}")
                        print(f"Target region: {top_candidate['target_region']}")
                        print(f"Usable length: {top_candidate['usable_length']} of {top_candidate['total_length']}")
                    
                    elif "target_sequence" in top_candidate:
                        print(f"Guide RNA: {top_candidate['target_sequence']}")
                        print(f"PAM: {top_candidate['pam_sequence']}")
                        print(f"Target region: {top_candidate['target_region']}")
                        print(f"Off-target risk: {top_candidate['off_target_risk']}")
            
            print("\nCONCLUSION:")
            print(report["conclusion"])
            print("="*80)
            
        else:
            print(f"Error: {message}")
            print("Please enter a valid DNA sequence containing only A, T, G, and C.")


if __name__ == "__main__":
    main()
