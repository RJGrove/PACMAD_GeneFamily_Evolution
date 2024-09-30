# -*- coding: utf-8 -*-
"""
Project: Deciphering Dynamics and Adaptations in PACMAD Gene Family Evolutionary Patterns
Script Function: Between running the BLASTN analysis and this script, it is necessary to run
                 the FASTA file containing only the conserved genes through MAFFT alignment
                 and then creating the most probable ancestor sequence through MEGA software.
                 After pasting the ancestor sequence into the FASTA containing the conserved
                 genes, align the entire sequence once more through MAFFT. Once completed,
                 paste the pathway to that aligned FASTA file into line 24. Results from
                 running this script will output a CSV file containing all identified
                 mutations labeled as well as several other columns of information.
                 
Author: RGrove
"""
############################################################################################
########################## Import necessary libraries and packages #########################
############################################################################################

from Bio import AlignIO
import pandas as pd
import matplotlib.pyplot as plt

gene_family_of_interest = "hom05m010550"

############################################################################################
############################################################################################
############################################################################################
""" Initialize species mapping dictionary and a placeholder for ancestor sequence """

alignment = AlignIO.read(f"C:\\Users\\rgrove4\\Documents\\RJGrove_SchnableLab_UCARE_PACMAD_GeneFamilies\\evolutionary_analysis\\{gene_family_of_interest}\\DNA_sequences_aligned\\aligned_ancestral_conserved_{gene_family_of_interest}.fasta", "fasta")
ancestor_id = "Most_Probable_Ancestor" 

species_mapping = {
    "Dexi": "Digitaria exilis",
    "Mis": "Miscanthus sinensis",
    "Seita": "Setaria italica",
    "Sevir": "Setaria viridis",
    "GWH": "Cenchrus purpureus",
    "Pahal": "Panicum hallii",
    "Ot": "Oropetium thomaeum",
    "Sspon": "Saccharum spontaneum",
    "Sobic": "Sorghum bicolor",
    "EJ": "Eragrostis curvula",
    "Zm": "Zea mays B73"
}

ancestor_sequence = None

############################################################################################
############################################################################################
############################################################################################
""" Assigns ancestor sequence from alignment; raises error if not found """

for record in alignment:
    if record.id == ancestor_id:
        ancestor_sequence = record
        break

if ancestor_sequence is None:
    raise ValueError("Ancestor sequence not found in alignment")

############################################################################################
############################################################################################
############################################################################################
""" List to input the mutation data into in order to save later on"""

mutation_data = []

############################################################################################
############################################################################################
############################################################################################
""" Identify species, initialize mutation tracking, and prepares to record positions """

for record in alignment:
    if record.id != ancestor_id:

        species_name = "Unknown Species"
        for prefix, full_name in species_mapping.items():
            if prefix in record.id:
                species_name = full_name
                break

        anc_seq = ancestor_sequence.seq
        current_seq = record.seq
        mutations = []
        
        deletion_positions = set()
        insertion_positions = set()

############################################################################################
############################################################################################
############################################################################################

        for i in range(min(len(anc_seq), len(current_seq))):
            anc_base = anc_seq[i]
            seq_base = current_seq[i]
            
            if anc_base != seq_base:
               
                # If the ancestor sequence has a gap and the current sequence has a base,
                # mark the position as an insertion.
                if anc_base == "-":
                    insertion_positions.add(i)
                
                # If the current sequence has a gap and the ancestor sequence has a base,
                # mark the position as a deletion.
                elif seq_base == "-":
                    deletion_positions.add(i)

############################################################################################
############################################################################################
############################################################################################

        for i in range(len(anc_seq)):
            
            # In case the sequences differ in length, assign a gap to those positions
            anc_base = anc_seq[i] if i < len(anc_seq) else "-"
            seq_base = current_seq[i] if i < len(current_seq) else "-"
            mutation_type = "No Mutation"
           
            if anc_base != seq_base:
               
                # If the ancestor has a gap and the current sequence has a base,
                # record as an insertion.
                if anc_base == "-" and seq_base != "-":
                    mutation_type = "Frameshift (Insertion)"
                
                # If the current sequence has a gap and the ancestor has a base,
                # record as a deletion.
                elif anc_base != "-" and seq_base == "-":
                    mutation_type = "Frameshift (Deletion)"
               
                # If both sequence have a base but they differ, record as a substitution.
                elif anc_base != "-" and seq_base != "-" and (anc_base != seq_base):
                    mutation_type = "Substitution"
                mutations.append([record.id, species_name, i + 1, anc_base, seq_base, True, mutation_type])
           
            else:
                mutations.append([record.id, species_name, i + 1, anc_base, seq_base, False, "No Mutation"])

############################################################################################
############################################################################################
############################################################################################

        for i in range(len(anc_seq)):
           
            # If the position is identified as an insertion or deletion, record "Indel"
            if i in insertion_positions or i in deletion_positions:
                mutation_type = "Indel"
                mutations.append([record.id, species_name, i + 1, anc_base, seq_base, True, mutation_type])

        mutation_data.extend(mutations)

############################################################################################
############################################################################################
############################################################################################
""" Save all additions to the dataframe into a csv file in current working directory """

df = pd.DataFrame(mutation_data, columns=["Gene ID", "Species Name", "Position", "Ancestor", "Current", "Mutation Present?", "Mutation Type"])
df.to_csv(f"{gene_family_of_interest}/datasets/mutations_analysis_{gene_family_of_interest}.csv", index=False)

############################################################################################
############################################################################################
############################################################################################
""" Plot the various mutations within each species as a bar graph """

plt.figure(figsize=(12, 8))
mutation_counts = df[df["Mutation Present?"] == True].groupby(["Species Name", "Mutation Type"]).size().unstack(fill_value=0)
mutation_counts.plot(kind='bar', stacked=True, figsize=(12, 8))
plt.title(f"Distribution of Mutation Types Across Species in {gene_family_of_interest}")
plt.xlabel("Species")
plt.ylabel("Count of Mutations")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Mutation Type")
plt.tight_layout()
plt.savefig(f'{gene_family_of_interest}/figures/{gene_family_of_interest}_mutation_distribution_bar_plot.png')
plt.show()

############################################################################################
############################################################################################
############################################################################################
""" Plot the percentages that mutations occur across all sequences in each species """

mutation_counts_per_species = df[df["Mutation Present?"] == True]["Species Name"].value_counts()
total_mutations = len(df[df["Mutation Present?"] == True])

plt.figure(figsize=(10, 7))
plt.pie(mutation_counts_per_species, labels=mutation_counts_per_species.index, autopct=lambda p: '{:.1f}%'.format(p), startangle=140, colors=plt.cm.Paired(range(len(mutation_counts_per_species))))
plt.title(f"Percentage of Mutations by Species in {gene_family_of_interest}")
plt.axis('equal')  
plt.savefig(f'{gene_family_of_interest}/figures/{gene_family_of_interest}_mutation_distribution_pie_chart.png')
plt.show()

############################################################################################
############################################################################################
############################################################################################