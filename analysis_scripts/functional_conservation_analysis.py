# -*- coding: utf-8 -*-
"""
Project: Deciphering Dynamics and Adaptations in PACMAD Gene Family Evolutionary Patterns
Script Function: 
    
Author: RGrove
"""
############################################################################################
########################## Import necessary libraries and packages #########################
############################################################################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

############################################################################################
############################################################################################
############################################################################################
""" Plot the distribution of percentage similarity and Shannon Entropy """

def plot_distributions(df):

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    
    axes[0].hist(df['percentage_similarity'], bins=20, color='blue', alpha=0.7)
    axes[0].set_title('Distribution of Percentage Similarity')
    axes[0].set_xlabel('Percentage Similarity')
    axes[0].set_ylabel('Frequency')
    
    axes[1].hist(df['Shannon_Entropy'], bins=20, color='green', alpha=0.7)
    axes[1].set_title('Distribution of Shannon Entropy')
    axes[1].set_xlabel('Shannon Entropy')
    axes[1].set_ylabel('Frequency')
    
    plt.tight_layout()
    plt.show()

############################################################################################
############################################################################################
############################################################################################
""" Plot the distribution of conserved gene families by their name """

def plot_conserved_gene_family_distribution(conserved_families):

    plt.figure(figsize=(10, 6))
    conserved_families['gene_family'].value_counts().plot(kind='bar', color='skyblue')
    plt.title('Conserved Genes in Homologous Gene Families for the PACMAD Clade')
    plt.xlabel('Gene Family')
    plt.ylabel('Number of Conserved Genes')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

############################################################################################
############################################################################################
############################################################################################
""" Identify conserved families and count species """

def find_conserved_gene_families(csv_file, similarity_threshold=20, entropy_threshold=1):
    
    df = pd.read_csv(csv_file)

    conserved_families = df[
        (df['percentage_similarity'] >= similarity_threshold) & 
        (df['Shannon_Entropy'] <= entropy_threshold)
    ]
    
    if conserved_families.empty:
        print(f"No conserved gene families found for thresholds: Similarity >= {similarity_threshold}, Entropy <= {entropy_threshold}")
        return None

    conserved_families['species_count'] = conserved_families.groupby('gene_family')['gene_id'].transform('count')
    conserved_families = conserved_families.sort_values(by='species_count', ascending=False)

    output_file = f'conserved_gene_families_similarity_{similarity_threshold}_entropy_{entropy_threshold}.csv'
    conserved_families.to_csv(output_file, index=False)
    
    print(f"Conserved gene families saved to {output_file}")

    return conserved_families[['gene_family', 'gene_family_type', 'gene_id', 'percentage_similarity', 'Shannon_Entropy', 'species_count']]

############################################################################################
############################################################################################
############################################################################################
""" Scatterplot of Shannon entropy (Variation) vs. Percent Similiarity (Conservation) """

def scatterplot(df):
    
    plot = sns.scatterplot(data = df, x = df['Shannon_Entropy'], y = df['percentage_similarity'], hue = df['gene_family'], marker = "x")
    plot.legend(fontsize = 7)
    
    plt.show()
    
############################################################################################
############################################################################################
############################################################################################

def main():
    csv_file = 'all_gene_family_sequence_analysis.csv'

    df = pd.read_csv(csv_file)

    plot_distributions(df)
    scatterplot(df)

    similarity_thresholds = [5, 10, 15, 20, 25]  
    entropy_thresholds = [1.0, 0.8, 0.5, 0.3]  

    for sim_thresh in similarity_thresholds:
        for ent_thresh in entropy_thresholds:
            print(f"Trying thresholds: Similarity >= {sim_thresh}, Entropy <= {ent_thresh}")
            conserved_families = find_conserved_gene_families(csv_file, sim_thresh, ent_thresh)
            
            if conserved_families is not None:
                print("Conserved Gene Families found!")
                print(conserved_families[['gene_family', 'gene_family_type', 'gene_id', 'percentage_similarity', 'Shannon_Entropy', 'species_count']])
                
                plot_conserved_gene_family_distribution(conserved_families)
                
                return

    print("No conserved gene families found with any of the tested thresholds.")

############################################################################################
############################################################################################
############################################################################################

if __name__ == "__main__":
    main()

############################################################################################
############################################################################################
############################################################################################
