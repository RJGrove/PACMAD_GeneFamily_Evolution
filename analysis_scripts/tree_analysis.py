# -*- coding: utf-8 -*-
"""
Project: Deciphering Dynamics and Adaptations in PACMAD Gene Family Evolutionary Patterns
Script Function: Once the variant analysis is performed, paste the pathway to that file in
                 line 19 of this script. Before running the code, make sure to input the
                 aligned FASTA file containing the ancestor and all genes sequences to MEGA
                 software to export a newick file containing the phylogenetic tree. Once
                 obtained, paste path to newick file in line 26 and then run the code. Results
                 will show a phylogenetic tree of all the genes and ancestor with highlights
                 to where mutations have been identified to occur across evolution.
                 
Author: RGrove
"""
############################################################################################
########################## Import necessary libraries and packages #########################
############################################################################################

import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle

gene_family_of_interest = "hom05m004593"

############################################################################################
############################################################################################
############################################################################################
""" Reads data, groups and analyzes mutations, and then loads the original tree """

mutation_df = pd.read_csv(f'{gene_family_of_interest}/datasets/mutations_analysis_{gene_family_of_interest}.csv')
mutation_counts = mutation_df.groupby('Species Name')['Mutation Present?'].sum()
conserved_mutations = mutation_df.groupby('Mutation Type')['Mutation Present?'].mean()
tree = Tree(f'{gene_family_of_interest}/phylogenetic_reconstructions/{gene_family_of_interest}_tree_reconstruction.nwk')

print("Leaf names (gene IDs) in the tree:")

############################################################################################
############################################################################################
############################################################################################
""" Highlights tree nodes based on gene IDs with mutations present """

for leaf in tree.iter_leaves():
    print(leaf.name)

for _, row in mutation_df.iterrows():
    if row['Mutation Present?']:

        print(f"Searching for gene ID: {row['Gene ID']}")
        
        nodes = tree.search_nodes(name=row['Gene ID'])
        
        if nodes:
            node = nodes[0]
            mutation_style = NodeStyle()
            mutation_style["bgcolor"] = "yellow"
            node.set_style(mutation_style)
        else:
            print(f"Gene ID not found in the tree: {row['Gene ID']}")

############################################################################################
############################################################################################
############################################################################################
""" Displays the phylogenetic tree with highlighted genes """

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True

tree.show(tree_style = ts)

############################################################################################
############################################################################################
############################################################################################