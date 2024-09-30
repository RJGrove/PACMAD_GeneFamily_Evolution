# -*- coding: utf-8 -*-
"""
Project: Deciphering Dynamics and Adaptations in PACMAD Gene Family Evolutionary Patterns
Script Function: 
    
Author: RGrove
"""
############################################################################################
########################## Import necessary libraries and packages #########################
############################################################################################

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import math

gene_family_of_interest = "hom05m010550" 

############################################################################################
############################################################################################
############################################################################################
""" Load aligned conserved gene sequences from a FASTA file """
    
def load_fasta(file_path):

    print(f"Loading aligned gene sequences from {file_path}")
    sequences = {}
   
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    
    print(f"Loaded {len(sequences)} sequences.")
    return sequences

############################################################################################
############################################################################################
############################################################################################
""" Analyze sequence variation compared to the Most Probable Ancestor and calculate the percentage of differences """

def analyze_sequence_variation(sequences, ancestor_seq):

    variation_results = {}
    ancestor_len = len(ancestor_seq)
    
    for gene_id, sequence in sequences.items():
        differences = sum(1 for a, b in zip(ancestor_seq, sequence) if a != b)
        percentage_difference = (differences / ancestor_len) * 100  # Calculate the percentage of differences
        percentage_similarity = (100 - percentage_difference)
        variation_results[gene_id] = {
            "differences_from_ancestor": differences,
            "percentage_difference": percentage_difference,
            "percentage_similarity": percentage_similarity
        }
    
    print(f"Analyzed variation for {len(sequences)} genes.")
    return variation_results

############################################################################################
############################################################################################
############################################################################################
""" Detect loss of function mutations (e.g., stop codons) in aligned sequences """

def detect_loss_of_function(sequences):
    
    loss_of_function_results = {}
    
    for gene_id, sequence in sequences.items():
        if '*' in sequence:
            loss_of_function_results[gene_id] = {"loss_of_function": "Stop codon detected"}
        else:
            loss_of_function_results[gene_id] = {"loss_of_function": "No loss detected"}
    
    print(f"Completed loss of function analysis for {len(sequences)} genes.")
    return loss_of_function_results

############################################################################################
############################################################################################
############################################################################################
""" Calculate the Shannon entropy for variability analysis """

def shannon_entropy(sequence):
    
    counts = Counter(sequence)
    total = len(sequence)
    entropy = -sum((count / total) * math.log2(count / total) for count in counts.values())
    return entropy

############################################################################################
############################################################################################
############################################################################################
""" Compile all results into a single DataFrame """

def compile_results(variation_results, loss_of_function_results, entropy_results):
    
    all_results = []
    
    for gene_id in variation_results.keys():
        result = {
            "gene_id": gene_id,
            "differences_from_ancestor": variation_results[gene_id]['differences_from_ancestor'],
            "percentage_difference": variation_results[gene_id]['percentage_difference'],
            "percentage_similarity": variation_results[gene_id]['percentage_similarity'],
            "loss_of_function": loss_of_function_results[gene_id]['loss_of_function'],
            "Shannon_Entropy": entropy_results.get(gene_id, np.nan)  # Use NaN if gene ID not found
        }
        all_results.append(result)

    return pd.DataFrame(all_results)

############################################################################################
############################################################################################
############################################################################################
""" Plot the percentage similarity of each gene to the ancestor sequence """

def plot_percentage_similarity(variation_results):
    
    gene_ids = list(variation_results.keys())
    similarities = [result['percentage_similarity'] for result in variation_results.values()]

    sorted_indices = np.argsort(similarities)[::-1]
    gene_ids_sorted = [gene_ids[i] for i in sorted_indices]
    similarities_sorted = [similarities[i] for i in sorted_indices]

    plt.figure(figsize=(12, 6))
    plt.bar(gene_ids_sorted, similarities_sorted, color='skyblue')
    plt.xlabel('Gene ID')
    plt.ylabel('Percentage Similarity to Ancestor')
    plt.title(f'Percentage Similarity of {gene_family_of_interest} Genes to Ancestor Sequence')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

############################################################################################
############################################################################################
############################################################################################

if __name__ == "__main__":
    aligned_fasta_file = f'C:\\Users\\rgrove4\\Documents\\RJGrove_SchnableLab_UCARE_PACMAD_GeneFamilies\\evolutionary_analysis\\{gene_family_of_interest}\\DNA_sequences_aligned\\aligned_ancestral_conserved_{gene_family_of_interest}.fasta'  
    aligned_sequences = load_fasta(aligned_fasta_file)
    
    ancestor_seq = "atggaggaggcggcggtgatggcggccggcgcccgtccctctggcgcgcgagccagcacggaggaggctgccagcgcccgcgccacgctcggcgcccgtccctccggcgcacgagccggcacggaggaggcggtcggcgcccacgcgcgcgcgcgctcggcgttcatccctccggcgcggcccacgtgcagcgaccgacacctaccagcgttctctctcttccaacgttcaaatattcgaaccgagcagcggttcgcgaggagccgaggaggtaaggaggaaaaggcactaggggcggcggcggcggcagaagcggcggcggcggccatggatttccacgccctctcgcgccgcgagctccaggcgctatgcaagcgcaacgccgtccgcgccaacatgtccaacgcggccatggcggacgcgctccggtccctcccctcggtgagcctatcctcgatcgatccccttcctttcctcctgagggtggacgggctcgacgagatcgggcgggcggtgcgggcgccgccgatgtcggcgatgaaatcggtggaggaggtgattatggatgaggagaagatcgacgggaacccgctcccccgcggcggccgtgcgcgttccaaggcaaggacggccgcgacggataagttggagcaggatgtaggagatcaagcccaagccaccttgcagggaatccaggggacggcggcaccacgggaggccttggcggcgcctctgaatgtggcaaaggaggtgaccagggaggagcaggggcatggttgcccgcttccccgcggtcgccgtgtaggggcgatgacaaggaaggctgcggcgcataagacaaagacggagggggatgaggaggaggaggac---gcggtgctggccccagcggacaccttgcagttccaggggagccgccagagaacggaggcaggggaggctggtgcggcccctgtggatgaggcagaggaggtagccacaggaaagagaaggacgacgaggtgtactagatcgaaggtgaagatggctttggatcagaaagaggaggtaccagcggcagcacagaaggagcagaaagtttcagacaagagctgcaatgatcctaaagaggatgaattggtagttgtggtggaggaagaggccacaaagccgcccgctcccccgcagtaaggaagatcagaatgtgacaaggaggactgcagtgcatgaggcagaggaggaggtgccggccccagccatcttgcggcggagccagcagaggacagtggcaccagagggcatggcacctgtggaggcggaagaggtagccacaggaaagagaaggacaaggaggtctgcgaggtcaaagacgaagatggctttggatcagaaagaggcagaagaggccacagaatgcaaggagcagaccaaggtgattcttctgatgtcgccactgggtctgtggtagtttcagaaaagagcatcgatgggtcctaaaacgcacgaagtggtagtaggcgtggtggtggaggaagatgtcatgaagccacaggaaggtcagaatgtgacaaggagacctccagtggatgagatggatgaggtaacgccagccatccctacagccatcctgaggcggagccaaaggacggcagcacgagaggctctggctcctgttgaggcaaaagaggtggccatggcaaagaaaatggcgagaaggtccacaagatctaaggt---gcaggtggctttggtcaagaaagaggaagcggcagcctcatttctctact-------------gcaaaacagaaggggcagaaagctgattcttctgttatgaccaatgggtctgctgtaatttataataagagcatcgattgtcctaaaaaggatgaagtagttgctgtagtggaggaagagaccacaaaaccacatgatggtgggaatgttaccaggaggatcacggcttatacacatgaaatggaggaggtgccggtcccagccacctcgcagcagagccaggggatggcagcactagaggctgcagtgcctgtgccacaggaaggtcaaaacgtgacgaggagagccgcaccacataagacacatgaggatttgcagcaagatcagagggaggcagaaaagggggctacaacaaagaggagaaggaggagaccaagaagatctaaagcagctgcagcagcacggaatgggcagaaagatttcgtggggacactgatgacgattgtttctgaatgtgccgctgaggatttgttgcttcttctgatgtaaccatagggtctgaagttgatcctaatgagcatgaagttgttcccaaggtggaggaggccgcaaaaccacaggaaggcgagaatgagacaaggaaggtcacagtgtataaggcacaggaggaggtgcagccaggccaggtgacagtggcaccagaggctatagcagctgtggaggcagaagaagtggccacggcaaagagaaagagaaggaggccaacaagatccaaggtgagggtggcggcatgtaaggaccagaaggccattgggtctcctgtattagccccagatcagagctgtgatgatcttagggaag---agcatcaacttattgcggcggtggaggaacatgtcacaaagccccaggaaggaattgtagaagaacaggatcaaagcagttctattcacaagtctgcttcctcagtaaaaatggaggatccaccaaccgtaagcacggtcactgagatggctccacctacaggctcagtccaagccagaagcaacatgtgttaattatatggcaacggctgctgatgaggagactgttaagatcaatcctgatactcaagacaagaaagactgtttcactctcaattctggggctggtcaattggattttcttgtcaatacacttaacagattttccaagcctatgcatgagttcaccatcaaggaggaggagaagaaggagggggagtgttggtggatgcttgtgtagactgccaacaaccatcaacagttaaacaaataa"
    variation_results = analyze_sequence_variation(aligned_sequences, ancestor_seq)
    loss_of_function_results = detect_loss_of_function(aligned_sequences)
    
    entropy_results = {}

    for gene_id, sequence in aligned_sequences.items():
        if gene_id != "Most_Probable_Ancestor":  
            entropy = shannon_entropy(sequence)
            entropy_results[gene_id] = entropy

    all_results_df = compile_results(variation_results, loss_of_function_results, entropy_results)
    all_results_df.to_csv(f'./{gene_family_of_interest}_sequence_analysis_results.csv', index=False)

    plot_percentage_similarity(variation_results)

    print(f"Analysis complete. Results saved to '{gene_family_of_interest}_sequence_analysis_results.csv'.")

############################################################################################
############################################################################################
############################################################################################