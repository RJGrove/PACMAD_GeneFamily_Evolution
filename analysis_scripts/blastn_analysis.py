"""
Project: Deciphering Dynamics and Adaptations in PACMAD Gene Family Evolutionary Patterns
Script Function: After running the BLASTN analysis on the gene family in question, this script
                 will perform several analyses in order to visualize the results. These analyses
                 consist of observing which genes are conserved over multiple species as well as 
                 an average count for the statisically values to better validate the findings.
                 Afterwards, the genes that are considered to be above threshold values will have
                 their DNA sequence data saved in a separate FASTA file.

Author: RGrove
"""
############################################################################################
########################## Import necessary libraries and packages #########################
############################################################################################

import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

gene_family_of_interest = "__________________" 

############################################################################################
############################################################################################
############################################################################################
""" Removing self-hits and applies an e-value threshold to BLASTN results """

df = pd.read_csv(f'{gene_family_of_interest}/blastn_analysis/{gene_family_of_interest}_blastn_output.csv')
filtered_df = df[df['query_id'] != df['subject_id']]

evalue_threshold = 2.9e-05 # Threshold was extracted from BLASTN to the banana genome; it's the maximum e-value
filtered_df = filtered_df[filtered_df['evalue'] <= evalue_threshold]
filtered_df.to_csv(f'{gene_family_of_interest}/blastn_analysis/filtered_{gene_family_of_interest}_blastn_output.csv', index = False)
 
############################################################################################
############################################################################################
############################################################################################
""" Counts conserved hits, filters by threshold, and calculates average metrics """

conservation_counts = filtered_df.groupby('query_id')['subject_id'].count()
threshold = 1  
conserved_across_species = conservation_counts[conservation_counts >= threshold]
print(conserved_across_species)

grouped = filtered_df.groupby('common_name').agg({
    'percent_identity': 'mean',
    'alignment_length': 'mean',
    'evalue': 'mean',
    'bit_score': 'mean'
})

############################################################################################
############################################################################################
############################################################################################
""" Plot bar graph of gene conservation and average metrics across gene family """

plt.figure(figsize = (12, 8))
plt.bar(conserved_across_species.index, conserved_across_species.values, color = 'gold')
plt.title(f'Genes Conserved Across {gene_family_of_interest}')
plt.xlabel('Gene ID')
plt.ylabel('No. of Comparisons')
plt.xticks(rotation = 90)
plt.grid(True)
plt.tight_layout()
plt.show()

def plot_bar_graph(column):
    plt.figure(figsize=(12, 8))
    plt.bar(grouped.index, grouped[column], color = 'skyblue' if column == 'percent_identity' else 'lightgreen' if column == 'alignment_length' else 'salmon')
    plt.title(f'Average {column} by Species for {gene_family_of_interest}')
    plt.xlabel('Species')
    plt.ylabel(column)
    plt.xticks(rotation=90)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plot_bar_graph('percent_identity')
plot_bar_graph('alignment_length')
plot_bar_graph('evalue')

############################################################################################
############################################################################################
############################################################################################
""" Creates list of conserved gene IDs and sets paths for input/output files """

conserved_genes = conserved_across_species.index.tolist()

original_fasta_path = f'{gene_family_of_interest}/DNA_sequences_raw/{gene_family_of_interest}.sequences.dna.fasta' # Change based on the original query fasta file
output_fasta_path = f'{gene_family_of_interest}/blastn_analysis/conserved_genes_sequences_{gene_family_of_interest}.fasta' # Change based on what the output fasta file will be named

############################################################################################
############################################################################################
############################################################################################
""" Extract sequences for only conserved genes and saves them into different FASTA file """

def extract_conserved_sequences(fasta_file, gene_ids, output_file):
    with open(fasta_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in gene_ids:
                SeqIO.write(record, outfile, "fasta")

extract_conserved_sequences(original_fasta_path, conserved_genes, output_fasta_path)
print(f"Conserved gene sequences saved to {output_fasta_path}")

############################################################################################
############################################################################################
############################################################################################
############################################################################################

