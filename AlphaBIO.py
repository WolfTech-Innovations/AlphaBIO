import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import random
import string
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import seaborn as sns
from collections import Counter
from Bio import SeqIO, Seq, SeqRecord, AlignIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.Blast import NCBIWWW, NCBIXML
from sklearn.metrics import accuracy_score, precision_recall_curve, auc
import os
import time
import json
import warnings
warnings.filterwarnings('ignore')

COMPANY_NAME = "WolfTech Innovations"
MODEL_NAME = "AlphaBIO-X"
VERSION = "3.0"

# ======================= ADVANCED DNA MODELING =======================

class DNAStructure:
    """Advanced DNA structure modeling with physical and chemical properties"""
    
    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self.length = len(sequence)
        self._compute_properties()
    
    def _compute_properties(self):
        """Calculate biophysical properties of the DNA sequence"""
        # Basic properties
        self.gc_content = self._calculate_gc_content()
        self.melting_temp = self._calculate_melting_temperature()
        self.molecular_weight = self._calculate_molecular_weight()
        
        # Secondary structure properties
        self.hairpins = self._find_hairpins()
        self.repeats = self._find_repeats()
        self.stability_profile = self._calculate_stability_profile()
        
        # Functional properties
        self.potential_orfs = self._find_open_reading_frames()
        self.restriction_sites = self._find_restriction_sites()
        self.modification_sites = self._find_modification_sites()
        
    def _calculate_gc_content(self):
        """Calculate GC content percentage"""
        gc_count = sum(1 for base in self.sequence if base in "GC")
        return (gc_count / self.length) * 100 if self.length > 0 else 0
    
    def _calculate_melting_temperature(self):
        """Calculate melting temperature using nearest-neighbor method"""
        # Simplified version of nearest-neighbor calculation
        if self.length < 14:
            # Wallace rule for short sequences
            return 2 * sum(1 for base in self.sequence if base in "AT") + 4 * sum(1 for base in self.sequence if base in "GC")
        else:
            # Advanced nearest-neighbor method (simplified)
            nearest_neighbors = {
                'AA': 9.1, 'AT': 8.6, 'TA': 6.0, 'CA': 5.8,
                'GT': 6.5, 'CT': 7.8, 'GA': 5.6, 'CG': 11.9,
                'GC': 11.1, 'GG': 11.0, 'TT': 9.1, 'TG': 5.8,
                'AC': 6.5, 'AG': 7.8, 'TC': 5.6, 'CC': 11.0
            }
            
            energy = sum(nearest_neighbors.get(self.sequence[i:i+2], 7.0) for i in range(len(self.sequence)-1))
            return 64.9 + 41 * (energy / (self.length - 1)) - (500 / self.length)
    
    def _calculate_molecular_weight(self):
        """Calculate molecular weight of the DNA sequence"""
        base_weights = {'A': 313.2, 'T': 304.2, 'G': 329.2, 'C': 289.2}
        return sum(base_weights.get(base, 0) for base in self.sequence) - 61.0  # Subtract water molecule
    
    def _find_hairpins(self):
        """Find potential hairpin structures in the DNA"""
        hairpins = []
        min_loop_size = 3
        min_stem_size = 3
        
        for i in range(self.length - (2*min_stem_size + min_loop_size)):
            for loop_size in range(min_loop_size, 10):
                for stem_size in range(min_stem_size, 8):
                    if i + stem_size + loop_size + stem_size > self.length:
                        continue
                    
                    stem1 = self.sequence[i:i+stem_size]
                    stem2 = self.sequence[i+stem_size+loop_size:i+stem_size+loop_size+stem_size]
                    stem2_rev_comp = self._reverse_complement(stem2)
                    
                    # Check if stems are complementary
                    matches = sum(1 for a, b in zip(stem1, stem2_rev_comp) if self._are_complementary(a, b))
                    if matches >= stem_size * 0.8:  # 80% match threshold
                        stability = matches / stem_size
                        position = i
                        hairpins.append({
                            'position': position,
                            'stem1': stem1,
                            'loop': self.sequence[i+stem_size:i+stem_size+loop_size],
                            'stem2': stem2,
                            'stability': stability
                        })
        
        return sorted(hairpins, key=lambda x: x['stability'], reverse=True)
    
    def _find_repeats(self):
        """Find direct and inverted repeats in the sequence"""
        repeats = {'direct': [], 'inverted': []}
        min_repeat_size = 5
        
        # Direct repeats
        for i in range(self.length - min_repeat_size):
            for j in range(i + min_repeat_size, self.length - min_repeat_size + 1):
                k = 0
                while i + k < self.length and j + k < self.length and self.sequence[i+k] == self.sequence[j+k]:
                    k += 1
                
                if k >= min_repeat_size:
                    repeats['direct'].append({
                        'seq': self.sequence[i:i+k],
                        'pos1': i,
                        'pos2': j,
                        'length': k
                    })
        
        # Inverted repeats
        for i in range(self.length - min_repeat_size):
            for j in range(i + min_repeat_size, self.length - min_repeat_size + 1):
                k = 0
                while (i + k < j and j + k < self.length and 
                       self._are_complementary(self.sequence[i+k], self.sequence[self.length-1-(j+k)])):
                    k += 1
                
                if k >= min_repeat_size:
                    repeats['inverted'].append({
                        'seq1': self.sequence[i:i+k],
                        'seq2': self.sequence[j:j+k],
                        'pos1': i,
                        'pos2': j,
                        'length': k
                    })
        
        return repeats
    
    def _calculate_stability_profile(self):
        """Calculate stability profile across the sequence"""
        if self.length < 2:
            return []
            
        # Free energy parameters for nearest neighbors (simplified)
        free_energy = {
            'AA': -1.0, 'AT': -0.9, 'TA': -0.6, 'CA': -1.3,
            'GT': -1.4, 'CT': -1.4, 'GA': -1.3, 'CG': -2.1,
            'GC': -2.3, 'GG': -1.8, 'TT': -1.0, 'TG': -1.3,
            'AC': -1.4, 'AG': -1.4, 'TC': -1.3, 'CC': -1.8
        }
        
        stability = []
        for i in range(self.length - 1):
            dinucleotide = self.sequence[i:i+2]
            energy = free_energy.get(dinucleotide, -1.0)
            stability.append(energy)
        
        return stability
    
    def _find_open_reading_frames(self):
        """Find potential open reading frames (ORFs)"""
        orfs = []
        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]
        min_orf_length = 30  # Minimum 10 amino acids
        
        for frame in range(3):
            i = frame
            while i < self.length - 2:
                if self.sequence[i:i+3] == start_codon:
                    start_pos = i
                    j = i + 3
                    while j < self.length - 2:
                        if self.sequence[j:j+3] in stop_codons:
                            if j - start_pos >= min_orf_length:
                                orfs.append({
                                    'start': start_pos,
                                    'stop': j + 2,
                                    'frame': frame,
                                    'length': j - start_pos + 3,
                                    'sequence': self.sequence[start_pos:j+3]
                                })
                            break
                        j += 3
                    i = j + 3 if j < self.length - 2 else self.length
                else:
                    i += 3
        
        return sorted(orfs, key=lambda x: x['length'], reverse=True)
    
    def _find_restriction_sites(self):
        """Find common restriction enzyme sites"""
        restriction_sites = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'XhoI': 'CTCGAG',
            'NotI': 'GCGGCCGC',
            'PstI': 'CTGCAG',
            'SmaI': 'CCCGGG',
            'SalI': 'GTCGAC'
        }
        
        sites = []
        for enzyme, site in restriction_sites.items():
            start = 0
            while True:
                start = self.sequence.find(site, start)
                if start == -1:
                    break
                sites.append({
                    'enzyme': enzyme,
                    'position': start,
                    'site': site
                })
                start += 1
        
        return sites
    
    def _find_modification_sites(self):
        """Find potential DNA modification sites"""
        # Methylation sites (CpG islands)
        cpg_sites = []
        for i in range(self.length - 1):
            if self.sequence[i:i+2] == "CG":
                cpg_sites.append(i)
        
        # Other modification sites
        other_sites = []
        motifs = {
            'GATC': 'Dam methylation',
            'GANTC': 'Hinf I',
            'GCGC': 'HhaI/HpaII methylation'
        }
        
        for motif, desc in motifs.items():
            if '*' in motif:
                # Handle ambiguous bases
                pass
            else:
                start = 0
                while True:
                    start = self.sequence.find(motif, start)
                    if start == -1:
                        break
                    other_sites.append({
                        'type': desc,
                        'position': start,
                        'motif': motif
                    })
                    start += 1
        
        return {'cpg': cpg_sites, 'other': other_sites}
    
    def _are_complementary(self, base1, base2):
        """Check if two bases are complementary"""
        complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return complements.get(base1, '') == base2
    
    def _reverse_complement(self, sequence):
        """Get the reverse complement of a DNA sequence"""
        complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complements.get(base, 'N') for base in reversed(sequence))
    
    def stability_score(self):
        """Calculate overall stability score (0-100)"""
        # Factors affecting stability
        factors = {
            'gc_content': 0.3,  # Weight for GC content
            'hairpins': 0.2,    # Weight for hairpin structures
            'repeats': 0.1,     # Weight for repeat regions
            'entropy': 0.2,     # Weight for sequence complexity
            'modification': 0.2  # Weight for modification potential
        }
        
        # GC content score (optimal around 50%)
        gc_score = 100 - abs(self.gc_content - 50) * 2
        
        # Hairpin score (fewer is better)
        hairpin_score = 100 - min(len(self.hairpins) * 10, 100)
        
        # Repeats score (fewer is better)
        repeat_count = len(self.repeats['direct']) + len(self.repeats['inverted'])
        repeat_score = 100 - min(repeat_count * 5, 100)
        
        # Sequence complexity/entropy score
        bases = Counter(self.sequence)
        entropy = -sum((count/self.length) * np.log2(count/self.length) for count in bases.values() if count > 0)
        max_entropy = np.log2(len(bases))
        entropy_score = (entropy / max_entropy) * 100 if max_entropy > 0 else 0
        
        # Modification potential score
        mod_count = len(self.modification_sites['cpg']) + len(self.modification_sites['other'])
        mod_density = mod_count / self.length if self.length > 0 else 0
        mod_score = 100 - min(mod_density * 1000, 100)  # Lower is better
        
        # Weighted average
        total_score = (
            gc_score * factors['gc_content'] +
            hairpin_score * factors['hairpins'] +
            repeat_score * factors['repeats'] +
            entropy_score * factors['entropy'] +
            mod_score * factors['modification']
        )
        
        return max(0, min(100, total_score))
    
    def binding_affinity(self, target_seq):
        """Calculate binding affinity to target sequence"""
        # Convert sequences to numpy arrays for vectorized operations
        seq1 = np.array(list(self.sequence))
        seq2 = np.array(list(self._reverse_complement(target_seq)))
        
        # Calculate lengths for alignment
        len1, len2 = len(seq1), len(seq2)
        min_len = min(len1, len2)
        
        # Perfect match would be all complementary bases
        # Calculate percentage match for different alignments
        best_score = 0
        best_offset = 0
        
        # Try different alignments
        for offset in range(-len1+1, len2):
            # Determine overlap region
            start1 = max(0, offset)
            end1 = min(len1, len2 + offset)
            start2 = max(0, -offset)
            end2 = min(len2, len1 - offset)
            
            overlap_len = end1 - start1
            
            if overlap_len <= 0:
                continue
            
            # Compare bases in overlap region
            slice1 = seq1[start1:end1]
            slice2 = seq2[start2:end2]
            
            # Calculate matches
            matches = {
                'AT_TA': np.sum(((slice1 == 'A') & (slice2 == 'T')) | ((slice1 == 'T') & (slice2 == 'A'))),
                'GC_CG': np.sum(((slice1 == 'G') & (slice2 == 'C')) | ((slice1 == 'C') & (slice2 == 'G')))
            }
            
            total_matches = matches['AT_TA'] + matches['GC_CG']
            match_percentage = (total_matches / overlap_len) * 100
            
            # Weight by overlap length as a percentage of the shorter sequence
            overlap_weight = overlap_len / min_len
            weighted_score = match_percentage * overlap_weight
            
            if weighted_score > best_score:
                best_score = weighted_score
                best_offset = offset
        
        # Calculate detailed binding statistics for best alignment
        if best_score > 0:
            start1 = max(0, best_offset)
            end1 = min(len1, len2 + best_offset)
            start2 = max(0, -best_offset)
            end2 = min(len2, len1 - best_offset)
            
            slice1 = seq1[start1:end1]
            slice2 = seq2[start2:end2]
            
            at_matches = np.sum(((slice1 == 'A') & (slice2 == 'T')) | ((slice1 == 'T') & (slice2 == 'A')))
            gc_matches = np.sum(((slice1 == 'G') & (slice2 == 'C')) | ((slice1 == 'C') & (slice2 == 'G')))
            
            return {
                'score': best_score,
                'offset': best_offset,
                'alignment_length': end1 - start1,
                'AT_matches': int(at_matches),
                'GC_matches': int(gc_matches),
                'total_matches': int(at_matches + gc_matches),
                'match_percentage': (at_matches + gc_matches) / (end1 - start1) * 100,
                'seq1_aligned': ''.join(seq1[start1:end1]),
                'seq2_aligned': ''.join(seq2[start2:end2])
            }
        else:
            return {'score': 0, 'offset': 0, 'alignment_length': 0}
    
    def effectiveness_against(self, target_dna):
        """Simulate effectiveness against target DNA"""
        binding = self.binding_affinity(target_dna)
        stability = self.stability_score()
        
        # Calculate overall effectiveness based on binding and stability
        binding_weight = 0.7  # Binding is more important
        stability_weight = 0.3  # Stability is also important
        
        # Effectiveness is a product of binding score and stability
        effectiveness = (binding['score'] * binding_weight) + (stability * stability_weight)
        
        # Bonus for longer binding regions
        min_binding_length = 15
        if binding['alignment_length'] > min_binding_length:
            length_bonus = min(20, (binding['alignment_length'] - min_binding_length) * 0.5)
            effectiveness += length_bonus
        
        # Bonus for GC-rich binding (stronger bonds)
        if binding['GC_matches'] > 0 and binding['alignment_length'] > 0:
            gc_percentage = binding['GC_matches'] / binding['alignment_length'] * 100
            gc_bonus = (gc_percentage - 50) * 0.1 if gc_percentage > 50 else 0
            effectiveness += gc_bonus
        
        # Cap at 100%
        return min(100, effectiveness)

# ======================= DNA SEQUENCE GENERATION =======================

class DNAGenerator:
    """Advanced DNA sequence generator with biological constraints"""
    
    def __init__(self):
        # Codon usage frequencies (human genome)
        self.codon_usage = {
            'TTT': 0.022, 'TTC': 0.018, 'TTA': 0.009, 'TTG': 0.013,
            'CTT': 0.014, 'CTC': 0.018, 'CTA': 0.008, 'CTG': 0.039,
            'ATT': 0.020, 'ATC': 0.020, 'ATA': 0.009, 'ATG': 0.022,
            'GTT': 0.013, 'GTC': 0.014, 'GTA': 0.008, 'GTG': 0.028,
            'TCT': 0.016, 'TCC': 0.017, 'TCA': 0.014, 'TCG': 0.004,
            'CCT': 0.018, 'CCC': 0.019, 'CCA': 0.018, 'CCG': 0.007,
            'ACT': 0.014, 'ACC': 0.018, 'ACA': 0.016, 'ACG': 0.006,
            'GCT': 0.018, 'GCC': 0.027, 'GCA': 0.016, 'GCG': 0.007,
            'TAT': 0.013, 'TAC': 0.015, 'TAA': 0.001, 'TAG': 0.001,
            'CAT': 0.011, 'CAC': 0.015, 'CAA': 0.013, 'CAG': 0.033,
            'AAT': 0.019, 'AAC': 0.018, 'AAA': 0.026, 'AAG': 0.031,
            'GAT': 0.021, 'GAC': 0.024, 'GAA': 0.030, 'GAG': 0.040,
            'TGT': 0.011, 'TGC': 0.012, 'TGA': 0.001, 'TGG': 0.012,
            'CGT': 0.005, 'CGC': 0.010, 'CGA': 0.006, 'CGG': 0.011,
            'AGT': 0.013, 'AGC': 0.019, 'AGA': 0.013, 'AGG': 0.012,
            'GGT': 0.011, 'GGC': 0.022, 'GGA': 0.017, 'GGG': 0.016
        }
        
        # Dinucleotide frequencies
        self.dinucleotide_freq = {
            'AA': 0.095, 'AT': 0.071, 'AG': 0.079, 'AC': 0.057,
            'TA': 0.065, 'TT': 0.095, 'TG': 0.082, 'TC': 0.063,
            'GA': 0.071, 'GT': 0.059, 'GG': 0.105, 'GC': 0.066,
            'CA': 0.071, 'CT': 0.075, 'CG': 0.018, 'CC': 0.086
        }
    
    def generate_random_dna(self, length=30, gc_content=0.5):
        """Generate random DNA sequence with specified GC content"""
        bases = ['G', 'C', 'A', 'T']
        probs = [gc_content/2, gc_content/2, (1-gc_content)/2, (1-gc_content)/2]
        
        return ''.join(np.random.choice(bases, size=length, p=probs))
    
    def generate_realistic_dna(self, length=30):
        """Generate biologically realistic DNA sequence using codon and dinucleotide frequencies"""
        if length < 2:
            return self.generate_random_dna(length)
        
        # Start with a random dinucleotide
        dinucleotides = list(self.dinucleotide_freq.keys())
        probabilities = list(self.dinucleotide_freq.values())
        
        # Normalize probabilities
        total = sum(probabilities)
        probabilities = [p/total for p in probabilities]
        
        # Generate first dinucleotide
        sequence = np.random.choice(dinucleotides, p=probabilities)
        
        # Markov chain approach to extend sequence
        while len(sequence) < length:
            last_base = sequence[-1]
            # Calculate conditional probabilities for next base
            next_base_probs = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
            
            total_prob = 0
            for dinuc, prob in self.dinucleotide_freq.items():
                if dinuc[0] == last_base:
                    next_base_probs[dinuc[1]] += prob
                    total_prob += prob
            
            # Normalize
            if total_prob > 0:
                for base in next_base_probs:
                    next_base_probs[base] /= total_prob
            
            # Select next base
            bases = list(next_base_probs.keys())
            probs = list(next_base_probs.values())
            next_base = np.random.choice(bases, p=probs)
            sequence += next_base
        
        return sequence
    
    def generate_targeted_dna(self, target_sequence, mutation_rate=0.1, indel_rate=0.05):
        """Generate DNA designed to target a specific sequence"""
        # Create complementary sequence as a starting point
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        complementary = ''.join(complement_map.get(base, 'N') for base in target_sequence)
        
        # Apply controlled mutations
        result = list(complementary)
        for i in range(len(result)):
            if random.random() < mutation_rate:
                original = result[i]
                # Avoid same base
                new_base = random.choice([b for b in "ATGC" if b != original])
                result[i] = new_base
        
        # Apply insertions and deletions (indels)
        if random.random() < indel_rate and len(result) > 10:
            # Deletion
            del_start = random.randint(0, len(result) - 3)
            del_length = random.randint(1, 2)  # Short deletions
            result = result[:del_start] + result[del_start + del_length:]
        
        if random.random() < indel_rate:
            # Insertion
            ins_pos = random.randint(0, len(result))
            ins_length = random.randint(1, 2)  # Short insertions
            insertion = [random.choice("ATGC") for _ in range(ins_length)]
            result = result[:ins_pos] + insertion + result[ins_pos:]
        
        return ''.join(result)
    
    def optimize_binding(self, target_sequence, population_size=50, generations=10):
        """Evolutionary approach to optimize binding to target sequence"""
        # Initial population
        population = []
        for _ in range(population_size):
            # Mix of complementary-based and random sequences
            if random.random() < 0.7:
                # Complementary-based with mutations
                individual = self.generate_targeted_dna(target_sequence, 
                                                       mutation_rate=random.uniform(0.05, 0.2))
            else:
                # Completely random
                individual = self.generate_realistic_dna(len(target_sequence))
            
            # Evaluate fitness
            dna_struct = DNAStructure(individual)
            fitness = dna_struct.effectiveness_against(target_sequence)
            population.append((individual, fitness))
        
        # Evolution process
        for generation in range(generations):
            # Sort by fitness
            population.sort(key=lambda x: x[1], reverse=True)
            
            # Keep top performers
            elite_size = max(2, population_size // 10)
            next_generation = population[:elite_size]
            
            # Create offspring for the rest of the population
            while len(next_generation) < population_size:
                # Selection (tournament selection)
                tournament_size = 3
                parent1 = max(random.sample(population, tournament_size), key=lambda x: x[1])[0]
                parent2 = max(random.sample(population, tournament_size), key=lambda x: x[1])[0]
                
                # Crossover
                crossover_point = random.randint(1, len(target_sequence) - 1)
                child = parent1[:crossover_point] + parent2[crossover_point:]
                
                # Mutation
                child_list = list(child)
                for i in range(len(child_list)):
                    if random.random() < 0.05:  # 5% mutation rate
                        child_list[i] = random.choice("ATGC")
                child = ''.join(child_list)
                
                # Evaluate fitness
                dna_struct = DNAStructure(child)
                fitness = dna_struct.effectiveness_against(target_sequence)
                
                next_generation.append((child, fitness))
            
            population = next_generation
        
        # Return best individual
        population.sort(key=lambda x: x[1], reverse=True)
        return population[0][0], population[0][1]
    
    def generate_disarming_dna(self, target_sequence, optimize_steps=3):
        """Generate multiple optimized sequences to disarm target DNA"""
        candidates = []
        
        # Multiple optimization runs with different parameters
        for _ in range(optimize_steps):
            pop_size = random.randint(30, 70)
            generations = random.randint(8, 15)
            sequence, fitness = self.optimize_binding(target_sequence, pop_size, generations)
            candidates.append((sequence, fitness))
        
        # Additional candidates using different approaches
        # 1. Perfect complementary as baseline
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        complementary = ''.join(complement_map.get(base, 'N') for base in target_sequence)
        compl_struct = DNAStructure(complementary)
        compl_fitness = compl_struct.effectiveness_against(target_sequence)
        candidates.append((complementary, compl_fitness))
        
        # 2. GC-enriched complementary
        gc_enriched = list(complementary)
        at_positions = [i for i, base in enumerate(gc_enriched) if base in "AT"]
        # Convert some A/T to G/C
        for pos in random.sample(at_positions, min(len(at_positions) // 3, 3)):
            gc_enriched[pos] = random.choice("GC")
        gc_enriched = ''.join(gc_enriched)
        gc_struct = DNAStructure(gc_enriched)
        gc_fitness = gc_struct.effectiveness_against(target_sequence)
        candidates.append((gc_enriched, gc_fitness))
        
        # Return best candidate
        candidates.sort(key=lambda x: x[1], reverse=True)
        return candidates[0][0], candidates[0][1], candidates

class AlphaBIO(nn.Module):
    """Advanced neural network for DNA analysis and cure generation"""
    
    def __init__(self, input_size=4, hidden_size=128, num_layers=3, dropout_rate=0.3):
        super(AlphaBIO, self).__init__()
        
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        
        # Embedding layer to convert DNA sequences to vector representations
        self.embedding = nn.Embedding(5, input_size)  # 4 nucleotides + padding
        
        # Bidirectional LSTM for sequence analysis
        self.lstm = nn.LSTM(
            input_size=input_size,
            hidden_size=hidden_size,
            num_layers=num_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout_rate if num_layers > 1 else 0
        )
        
        # Attention mechanism
        self.attention = nn.Sequential(
            nn.Linear(hidden_size * 2, hidden_size),
            nn.Tanh(),
            nn.Linear(hidden_size, 1),
            nn.Softmax(dim=1)
        )
        
        # Fully connected layers for classification and regression tasks
        self.fc_classify = nn.Sequential(
            nn.Linear(hidden_size * 2, hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_size, hidden_size // 2),
            nn.ReLU(),
            nn.Linear(hidden_size // 2, 2)  # Binary classification (effective/not effective)
        )
        
        # Regression head for binding affinity prediction
        self.fc_affinity = nn.Sequential(
            nn.Linear(hidden_size * 2, hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_size, 1),
            nn.Sigmoid()  # Output between 0-1 (will be scaled to percentage)
        )
        
        # Sequence generation head
        self.fc_generate = nn.Sequential(
            nn.Linear(hidden_size * 2, hidden_size * 4),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_size * 4, hidden_size * 8),
            nn.ReLU(),
            nn.Linear(hidden_size * 8, 4)  # Logits for 4 nucleotides
        )
        
        # Initialize weights
        self._init_weights()
    
    def _init_weights(self):
        """Initialize weights for better convergence"""
        for name, param in self.named_parameters():
            if 'weight' in name:
                nn.init.xavier_normal_(param)
            elif 'bias' in name:
                nn.init.constant_(param, 0.0)
    
    def _encode_sequence(self, sequence):
        """Convert DNA sequence string to tensor"""
        # Map nucleotides to indices
        nucleotide_map = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'N': 4}
        return torch.tensor([nucleotide_map.get(nuc, 4) for nuc in sequence])
    
    def _apply_attention(self, lstm_output):
        """Apply attention mechanism to focus on important regions"""
        attention_weights = self.attention(lstm_output)
        context_vector = torch.sum(attention_weights * lstm_output, dim=1)
        return context_vector, attention_weights
    
    def forward(self, x, task='classify'):
        """
        Forward pass through the network
        
        Args:
            x: Input DNA sequence (can be string or tensor)
            task: One of ['classify', 'affinity', 'generate']
            
        Returns:
            Task-specific output
        """
        # Convert string input if needed
        if isinstance(x, str):
            x = self._encode_sequence(x).unsqueeze(0)  # Add batch dimension
        
        # Handle batch or single input
        if x.dim() == 1:
            x = x.unsqueeze(0)  # Add batch dimension
            
        # Embed sequence
        embedded = self.embedding(x)
        
        # Process through LSTM
        lstm_output, (hidden, cell) = self.lstm(embedded)
        
        # Apply attention
        context, attention_weights = self._apply_attention(lstm_output)
        
        # Branch based on task
        if task == 'classify':
            return self.fc_classify(context)
        elif task == 'affinity':
            return self.fc_affinity(context) * 100  # Scale to percentage
        elif task == 'generate':
            # For generation, we need logits for each position
            return self.fc_generate(lstm_output)
        else:
            raise ValueError(f"Unknown task: {task}")
    
    def predict_effectiveness(self, sequence, threshold=0.5):
        """Predict if a DNA sequence will be effective"""
        with torch.no_grad():
            logits = self.forward(sequence, task='classify')
            probs = torch.softmax(logits, dim=1)
            effective_prob = probs[0, 1].item()
            prediction = effective_prob > threshold
            return {
                'effective': prediction,
                'probability': effective_prob,
                'confidence': max(effective_prob, 1 - effective_prob)
            }
    
    def predict_binding_affinity(self, sequence, target_sequence):
        """Predict binding affinity between two sequences"""
        # Concatenate sequences with special separator
        combined = sequence + "N" + target_sequence
        
        with torch.no_grad():
            affinity = self.forward(combined, task='affinity').item()
            return affinity
    
    def generate_cure_sequence(self, target_sequence, length=30, temperature=0.8):
        """Generate a therapeutic DNA sequence to counter target sequence"""
        # Start with special "start" token
        current_sequence = "N"
        
        # Generate sequence one nucleotide at a time
        for _ in range(length):
            with torch.no_grad():
                # Get predictions for next nucleotide
                combined = current_sequence + "N" + target_sequence
                logits = self.forward(combined, task='generate')
                
                # Get probabilities for the last position
                last_pos_logits = logits[0, -1, :]
                
                # Apply temperature to control randomness
                scaled_logits = last_pos_logits / temperature
                probs = torch.softmax(scaled_logits, dim=0)
                
                # Sample from distribution
                next_nucleotide_idx = torch.multinomial(probs, 1).item()
                next_nucleotide = "ATGC"[next_nucleotide_idx]
                
                # Add to sequence
                current_sequence += next_nucleotide
        
        # Remove the start token
        return current_sequence[1:]
    
    def optimize_cure(self, target_sequence, population_size=10, generations=5):
        """Use genetic algorithm to optimize cure sequence"""
        # Generate initial population
        population = []
        for _ in range(population_size):
            sequence = self.generate_cure_sequence(target_sequence)
            affinity = self.predict_binding_affinity(sequence, target_sequence)
            effectiveness = self.predict_effectiveness(sequence)
            fitness = affinity * (0.5 + 0.5 * effectiveness['probability'])
            population.append((sequence, fitness))
        
        # Evolution process
        for generation in range(generations):
            # Sort by fitness
            population.sort(key=lambda x: x[1], reverse=True)
            
            # Keep top performers
            elite_size = max(1, population_size // 5)
            next_generation = population[:elite_size]
            
            # Create offspring for the rest of the population
            while len(next_generation) < population_size:
                # Selection (tournament selection)
                parent1 = max(random.sample(population, 3), key=lambda x: x[1])[0]
                parent2 = max(random.sample(population, 3), key=lambda x: x[1])[0]
                
                # Crossover
                crossover_point = random.randint(1, len(parent1) - 1)
                child = parent1[:crossover_point] + parent2[crossover_point:]
                
                # Mutation
                child_list = list(child)
                for i in range(len(child_list)):
                    if random.random() < 0.1:  # 10% mutation rate
                        child_list[i] = random.choice("ATGC")
                child = ''.join(child_list)
                
                # Evaluate fitness
                affinity = self.predict_binding_affinity(child, target_sequence)
                effectiveness = self.predict_effectiveness(child)
                fitness = affinity * (0.5 + 0.5 * effectiveness['probability'])
                
                next_generation.append((child, fitness))
            
            population = next_generation
        
        # Return best solution
        population.sort(key=lambda x: x[1], reverse=True)
        best_solution = population[0][0]
        
        return {
            'sequence': best_solution,
            'binding_affinity': self.predict_binding_affinity(best_solution, target_sequence),
            'effectiveness': self.predict_effectiveness(best_solution),
            'similarity': self._sequence_similarity(best_solution, target_sequence)
        }
    
    def _sequence_similarity(self, seq1, seq2):
        """Calculate sequence similarity between two sequences of potentially different lengths"""
        # Use the shorter sequence to calculate maximum possible matches
        min_len = min(len(seq1), len(seq2))
        matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
        return (matches / min_len) * 100 if min_len > 0 else 0

# ======================= SEQUENCE ANALYSIS AND CURE GENERATION =======================

class DNACureGenerator:
    """Class to coordinate DNA analysis and cure generation"""
    
    def __init__(self, model_path=None):
        # Initialize the model
        self.model = AlphaBIO()
        if model_path and os.path.exists(model_path):
            self.model.load_state_dict(torch.load(model_path))
        self.model.eval()  # Set to evaluation mode
        
        # DNA structure analyzer
        self.dna_analyzer = DNAStructure
        
        # Sequence generator
        self.sequence_generator = DNAGenerator()
        
        # Optimization settings
        self.optimization_rounds = 3
        self.advanced_mode = True
    
    def analyze_target(self, target_sequence):
        """Analyze a target DNA sequence to understand its properties"""
        # Basic structure analysis
        structure = self.dna_analyzer(target_sequence)
        stability = structure.stability_score()
        
        # Identify key features
        features = {
            'gc_content': structure.gc_content,
            'length': structure.length,
            'melting_temp': structure.melting_temp,
            'hairpins': len(structure.hairpins),
            'repeats': len(structure.repeats['direct']) + len(structure.repeats['inverted']),
            'orfs': len(structure.potential_orfs),
            'stability': stability,
            'modification_sites': len(structure.modification_sites['cpg']) + 
                                 len(structure.modification_sites['other'])
        }
        
        # Identify vulnerability points
        vulnerabilities = []
        
        # Check for regions with low stability
        stability_profile = structure.stability_profile
        for i in range(len(stability_profile) - 5):
            window = stability_profile[i:i+5]
            if sum(window) / 5 > -1.0:  # High values indicate instability
                vulnerabilities.append({
                    'type': 'low_stability',
                    'position': i,
                    'length': 5,
                    'score': sum(window) / 5
                })
        
        # Check for hairpin regions
        for hairpin in structure.hairpins[:3]:  # Top 3 most stable hairpins
            vulnerabilities.append({
                'type': 'hairpin',
                'position': hairpin['position'],
                'length': len(hairpin['stem1']) * 2 + len(hairpin['loop']),
                'score': hairpin['stability']
            })
        
        # Neural network analysis
        if self.advanced_mode:
            with torch.no_grad():
                effectiveness = self.model.predict_effectiveness(target_sequence)
                features['model_score'] = effectiveness['probability']
        
        return {
            'features': features,
            'vulnerabilities': sorted(vulnerabilities, key=lambda x: x['score'], reverse=True),
            'structure': structure
        }
    
    def generate_cure_candidates(self, target_sequence, num_candidates=5):
        """Generate multiple candidate cures for a target sequence"""
        candidates = []
        
        # Method 1: Using neural network (advanced mode)
        if self.advanced_mode:
            for _ in range(max(1, num_candidates // 3)):
                cure = self.model.optimize_cure(target_sequence)
                candidates.append({
                    'sequence': cure['sequence'],
                    'method': 'neural_network',
                    'binding_affinity': cure['binding_affinity'],
                    'effectiveness': cure['effectiveness']['probability']
                })
        
        # Method 2: Using evolutionary approach
        for _ in range(max(1, num_candidates // 3)):
            sequence, fitness, _ = self.sequence_generator.generate_disarming_dna(target_sequence)
            
            # Evaluate with neural network if available
            if self.advanced_mode:
                effectiveness = self.model.predict_effectiveness(sequence)
                binding = self.model.predict_binding_affinity(sequence, target_sequence)
            else:
                structure = self.dna_analyzer(sequence)
                effectiveness = {'probability': structure.effectiveness_against(target_sequence) / 100}
                binding_info = structure.binding_affinity(target_sequence)
                binding = binding_info['score']
            
            candidates.append({
                'sequence': sequence,
                'method': 'evolutionary',
                'binding_affinity': binding,
                'effectiveness': effectiveness['probability']
            })
        
        # Method 3: Direct complementary with optimizations
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        complementary = ''.join(complement_map.get(base, 'N') for base in target_sequence)
        
        # Evaluate with neural network if available
        if self.advanced_mode:
            effectiveness = self.model.predict_effectiveness(complementary)
            binding = self.model.predict_binding_affinity(complementary, target_sequence)
        else:
            structure = self.dna_analyzer(complementary)
            effectiveness = {'probability': structure.effectiveness_against(target_sequence) / 100}
            binding_info = structure.binding_affinity(target_sequence)
            binding = binding_info['score']
        
        candidates.append({
            'sequence': complementary,
            'method': 'complementary',
            'binding_affinity': binding,
            'effectiveness': effectiveness['probability']
        })
        
        # Sort candidates by predicted effectiveness
        candidates.sort(key=lambda x: x['effectiveness'], reverse=True)
        
        return candidates[:num_candidates]
    
    def evaluate_cure(self, cure_sequence, target_sequence):
        """Comprehensive evaluation of a cure against target sequence"""
        # Structure analysis
        cure_structure = self.dna_analyzer(cure_sequence)
        target_structure = self.dna_analyzer(target_sequence)
        
        # Basic binding analysis
        binding_info = cure_structure.binding_affinity(target_sequence)
        
        # Advanced analysis with neural network
        if self.advanced_mode:
            effectiveness = self.model.predict_effectiveness(cure_sequence)
            binding_score = self.model.predict_binding_affinity(cure_sequence, target_sequence)
        else:
            effectiveness = {
                'effective': binding_info['score'] > 70,
                'probability': binding_info['score'] / 100
            }
            binding_score = binding_info['score']
        
        # Calculate stability under different conditions
        stability = {
            'normal': cure_structure.stability_score(),
            'high_temp': self._simulate_stability(cure_sequence, temperature=45),
            'low_ph': self._simulate_stability(cure_sequence, ph=5.0)
        }
        
        # Estimated effectiveness percentage
        overall_score = (
            binding_score * 0.6 +  # 60% weight to binding
            effectiveness['probability'] * 100 * 0.3 +  # 30% weight to model prediction
            stability['normal'] * 0.1  # 10% weight to stability
        )
        
        return {
            'sequence': cure_sequence,
            'target': target_sequence,
            'length': len(cure_sequence),
            'gc_content': cure_structure.gc_content,
            'binding_data': binding_info,
            'binding_score': binding_score,
            'effectiveness': effectiveness,
            'stability': stability,
            'overall_score': min(100, overall_score),
            'visualization': self._visualize_binding(cure_sequence, target_sequence, binding_info)
        }
    
    def _simulate_stability(self, sequence, temperature=37, ph=7.0):
        """Simulate stability under different conditions"""
        structure = self.dna_analyzer(sequence)
        base_stability = structure.stability_score()
        
        # Temperature effect (stability decreases with higher temperature)
        if temperature > 37:
            temp_factor = 1.0 - min(1.0, (temperature - 37) / 50)
        else:
            temp_factor = 1.0 + min(0.2, (37 - temperature) / 50)
        
        # pH effect (stability is best around neutral pH)
        ph_deviation = abs(ph - 7.0)
        ph_factor = 1.0 - min(0.5, ph_deviation / 4.0)
        
        return base_stability * temp_factor * ph_factor
    
    def _visualize_binding(self, cure_sequence, target_sequence, binding_info):
        """Create a text-based visualization of the binding"""
        if 'seq1_aligned' not in binding_info or 'seq2_aligned' not in binding_info:
            return "Binding visualization not available"
        
        seq1 = binding_info['seq1_aligned']
        seq2 = binding_info['seq2_aligned']
        
        # Create binding representation
        matches = ""
        for a, b in zip(seq1, seq2):
            if (a == 'A' and b == 'T') or (a == 'T' and b == 'A') or (a == 'G' and b == 'C') or (a == 'C' and b == 'G'):
                matches += "|"
            else:
                matches += " "
        
        return f"Cure:   5'-{seq1}-3'\n       {matches}\nTarget: 3'-{seq2}-5'"
    
    def optimize_cure(self, initial_cure, target_sequence, iterations=5):
        """Optimize an initial cure to improve its effectiveness"""
        best_cure = initial_cure
        best_score = self.evaluate_cure(best_cure, target_sequence)['overall_score']
        
        for iteration in range(iterations):
            # Generate variations
            variations = []
            
            # Variation 1: Point mutations
            mutated = list(best_cure)
            mutation_count = random.randint(1, max(1, len(best_cure) // 10))
            for _ in range(mutation_count):
                pos = random.randint(0, len(mutated) - 1)
                original = mutated[pos]
                mutated[pos] = random.choice([b for b in "ATGC" if b != original])
            variations.append(''.join(mutated))
            
            # Variation 2: Local optimization at binding regions
            cure_structure = self.dna_analyzer(best_cure)
            binding_info = cure_structure.binding_affinity(target_sequence)
            
            if 'offset' in binding_info and binding_info['alignment_length'] > 0:
                # Target the binding region
                start = max(0, binding_info['offset']) if binding_info['offset'] >= 0 else 0
                end = min(start + binding_info['alignment_length'], len(best_cure))
                
                # Leave the binding region intact but optimize surrounding regions
                if start > 0:
                    prefix = self.sequence_generator.generate_realistic_dna(start)
                else:
                    prefix = ""
                
                if end < len(best_cure):
                    suffix = self.sequence_generator.generate_realistic_dna(len(best_cure) - end)
                else:
                    suffix = ""
                
                optimized = prefix + best_cure[start:end] + suffix
                variations.append(optimized)
            
            # Variation 3: Use the model to generate a new sequence (if available)
            if self.advanced_mode:
                model_generated = self.model.generate_cure_sequence(target_sequence, len(best_cure))
                variations.append(model_generated)
            
            # Evaluate all variations
            for variant in variations:
                score = self.evaluate_cure(variant, target_sequence)['overall_score']
                if score > best_score:
                    best_cure = variant
                    best_score = score
        
        # Final evaluation
        return self.evaluate_cure(best_cure, target_sequence)
    
    def get_report(self, target_sequence, cure_sequence=None):
        """Generate comprehensive report on target and potential cures"""
        # Analyze target
        target_analysis = self.analyze_target(target_sequence)
        
        # Generate or evaluate cure
        if cure_sequence:
            cure_evaluation = self.evaluate_cure(cure_sequence, target_sequence)
            cure_candidates = [cure_evaluation]
        else:
            # Generate candidates and evaluate best one
            candidates = self.generate_cure_candidates(target_sequence, num_candidates=3)
            cure_candidates = []
            for candidate in candidates:
                evaluation = self.evaluate_cure(candidate['sequence'], target_sequence)
                cure_candidates.append(evaluation)
        
        # Sort candidates by score
        cure_candidates.sort(key=lambda x: x['overall_score'], reverse=True)
        
        # Create report
        report = {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'model_info': {
                'name': MODEL_NAME,
                'version': VERSION,
                'company': COMPANY_NAME,
                'advanced_mode': self.advanced_mode
            },
            'target_analysis': target_analysis,
            'cure_candidates': cure_candidates,
            'recommendations': self._generate_recommendations(target_analysis, cure_candidates)
        }
        
        return report
    
    def _generate_recommendations(self, target_analysis, cure_candidates):
        """Generate recommendations based on analysis"""
        recommendations = []
        
        # Best cure recommendation
        if cure_candidates:
            best_cure = cure_candidates[0]
            if best_cure['overall_score'] >= 85:
                confidence = "high"
            elif best_cure['overall_score'] >= 70:
                confidence = "medium"
            else:
                confidence = "low"
                
            recommendations.append({
                'type': 'primary_cure',
                'confidence': confidence,
                'description': f"The primary recommended cure has an effectiveness score of {best_cure['overall_score']:.1f}%.",
                'sequence': best_cure['sequence']
            })
        
        # Target vulnerability recommendation
        vulnerabilities = target_analysis['vulnerabilities']
        if vulnerabilities:
            vulnerable_regions = []
            for v in vulnerabilities[:2]:  # Top 2 vulnerabilities
                start = v['position']
                end = v['position'] + v['length']
                vulnerable_regions.append(f"{start}-{end}")
                
            recommendations.append({
                'type': 'vulnerability',
                'description': f"Target DNA has key vulnerabilities in regions {', '.join(vulnerable_regions)}.",
                'vulnerable_regions': vulnerable_regions
            })
        
        # Stability recommendation
        if cure_candidates and cure_candidates[0]['stability']['normal'] < 70:
            recommendations.append({
                'type': 'stability_warning',
                'description': "The generated cure has below-optimal stability. Consider testing in controlled conditions."
            })
        
        return recommendations
        
def main():
    # Define a fallback dictionary with default values to avoid KeyError
    binding = {
        'GC_matches': 0,  # No GC matches
        'alignment_length': 0,  # No alignment length
        'binding_score': 0,  # Default binding score
    }
    # Get the sequence from user input
    sequence = input("Please enter a DNA sequence: ").strip()

    # Ensure the sequence is not empty or invalid (basic validation)
    if not sequence or not all(base in 'ATGC' for base in sequence.upper()):
        print("Invalid sequence. Please enter a valid DNA sequence consisting of A, T, G, and C.")
        return

    # Initialize the DNACureGenerator with the given model path (if any)
    cure_generator = DNACureGenerator()

    # Analyze the target sequence
    target_analysis = cure_generator.analyze_target(sequence)

    # Print the target sequence analysis
    print("\nTarget Sequence Analysis:")
    print(f"GC Content: {target_analysis['features']['gc_content']}%")
    print(f"Length: {target_analysis['features']['length']} base pairs")
    print(f"Melting Temperature: {target_analysis['features']['melting_temp']}°C")
    print(f"Stability Score: {target_analysis['features']['stability']}")
    print(f"Hairpins Found: {target_analysis['features']['hairpins']}")
    print(f"Repeats Found: {target_analysis['features']['repeats']} total")
    print(f"Open Reading Frames (ORFs): {target_analysis['features']['orfs']}")
    print(f"Modification Sites: {target_analysis['features']['modification_sites']}")
    
    # Print vulnerabilities found
    print("\nTarget Sequence Vulnerabilities:")
    for vulnerability in target_analysis['vulnerabilities']:
        print(f"Type: {vulnerability['type']}, Position: {vulnerability['position']}, Score: {vulnerability['score']}")

    # Generate cure candidates based on the target sequence
    print("\nGenerating Cure Candidates...")
    cure_candidates = cure_generator.generate_cure_candidates(sequence)

    # Evaluate and display the top cure candidates
    print("\nCure Candidate Evaluation:")
    for i, candidate in enumerate(cure_candidates):
        cure_evaluation = cure_generator.evaluate_cure(candidate['sequence'], sequence)
        print(f"\nCandidate {i + 1}:")
        print(f"Sequence: {candidate['sequence']}")
        print(f"Effectiveness: {cure_evaluation['effectiveness']['probability']*100:.2f}%")
        print(f"Binding Score: {cure_evaluation['binding_score']}")
        print(f"Stability: {cure_evaluation['stability']['normal']}/100")
        print(f"Overall Score: {cure_evaluation['overall_score']:.2f}%")
        
    # Optionally, select the best candidate for further optimization
    print("\nOptimizing Best Cure Candidate...")
    best_candidate = cure_candidates[0]
    optimized_cure = cure_generator.optimize_cure(best_candidate['sequence'], sequence)

    print(f"\nOptimized Cure: {optimized_cure['sequence']}")
    print(f"Optimized Effectiveness: {optimized_cure['effectiveness']['probability']*100:.2f}%")
    print(f"Optimized Binding Score: {optimized_cure['binding_score']}")
    print(f"Optimized Stability: {optimized_cure['stability']['normal']}/100")
    print(f"Optimized Overall Score: {optimized_cure['overall_score']:.2f}%")

if __name__ == "__main__":
    main()


