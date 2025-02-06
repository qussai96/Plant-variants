#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import gc
from Bio.Seq import Seq
import re
import os
from Bio import pairwise2
from collections import defaultdict
import math
import time
import sys

 
# ## Data can be downloaded from Ensembl Plant at https://plants.ensembl.org/info/data/ftp/index.html

# Input file with the mutation information in gvf or interchangeable form 
gvf_file = r'/proj/q.abbas/Plant_mutations/oryza_sativa_incl_consequences.gvf'

# Input fore the introns and CDS in gtf or interchangeable form. the 
gtf_file = r'/proj/q.abbas/Plant_mutations/Oryza_sativa.IRGSP-1.0.60.gtf'

# Input for the protein sequence data in FASTA format 
pep_file = r'/proj/q.abbas/Plant_mutations/Oryza_sativa.IRGSP-1.0.pep.all.fa'

# Input for the cds sequence data in FASTA format
cds_file = r'/proj/q.abbas/Plant_mutations/Oryza_sativa.IRGSP-1.0.cds.all.fa'


# Helper function to extract information from the 'Suppl' column
def extract_id(suppl_value, pattern):
    match = re.search(pattern, suppl_value)
    if match:
        return match.group(1)
    return None


def gtf_to_df_with_genes(gtf_file):
    column_names = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']
    df = pd.read_csv(gtf_file, sep='\t', index_col=False, names=column_names, dtype={'Start': int, 'End': int}, comment='#')
    df = df[df['Type'].isin(['gene', 'CDS'])]
    df['Genes'] = df.apply(lambda row: extract_id(row['Suppl'], r'transcript_id "([^"]+)"'), axis=1)
    return df


def add_introns(gtf_file):
    df = gtf_to_df_with_genes(gtf_file)
    output_dfs = []
    
    for _, sub_df in df.groupby('Genes'):
        df_cdss = sub_df[sub_df['Type'] == 'CDS'].sort_values('Start')
        if len(df_cdss) > 1:
            end_values = df_cdss['End'].values
            start_values = df_cdss['Start'].values
            intron_mask = start_values[1:] - end_values[:-1] > 1
            intron_regions = list(zip(end_values[:-1][intron_mask] + 1, start_values[1:][intron_mask] - 1))
        else:
            intron_regions = []
        
        rows = []
        for _, row in sub_df.iterrows():
            if row['Type'] == 'CDS':
                current_cds = (row['Start'], row['End'])
                for intron_region in intron_regions:
                    if current_cds[1] + 1 == intron_region[0]:
                        new_row = row.copy()
                        new_row['Type'] = 'intron'
                        new_row['Start'] = intron_region[0]
                        new_row['End'] = intron_region[1]
                        rows.append(row)
                        rows.append(new_row)
            else:
                rows.append(row)
        
        if len(df_cdss) > 0:
            last_cds = df_cdss.iloc[-1]
            if last_cds['End'] + 1 not in [region[0] for region in intron_regions]:
                rows.append(last_cds)
        
        output_dfs.append(pd.DataFrame(rows))
    
    output_df = pd.concat(output_dfs, ignore_index=True)
    output_df = output_df[['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    
    return output_df


# Helper function to clean the variant effect column. leave only the  variant effect without any additional information  
def clean_variant_effect(x):
    if isinstance(x, str):
        # Split the string by commas to handle multiple variant effects
        effects = x.split(',')
        
        # Clean each individual effect by removing everything after the first space
        cleaned_effects = [re.sub(r'\s+.*', '', effect).strip() for effect in effects]
        
        # Join the cleaned effects back with commas
        x = ','.join(cleaned_effects)
        
    return x


# Helper function to remove all variant peptide without a matching Transcript_ID from the variant_peptide column 
def filter_variant_peptide_by_transcript(df, variant_peptide_col='variant_peptide', transcript_col='Transcript_ID'):
    def keep_matching_variant_peptides(variant_peptides, transcript_id):
        # If the variant_peptides is NaN or empty, return an empty string
        if not variant_peptides:
            return ''
        
        # Split the variant_peptides by commas, then filter by whether the transcript_id is in each part
        filtered_peptides = [peptide for peptide in variant_peptides.split(',') if transcript_id in peptide]
        
        # Rejoin the filtered peptides into a single string or return an empty string if none match
        return ','.join(filtered_peptides) if filtered_peptides else ''

    # Apply the filtering to the DataFrame
    df[variant_peptide_col] = df.apply(
        lambda row: keep_matching_variant_peptides(row[variant_peptide_col], row[transcript_col]), axis=1
    )

    return df


# Helper function to explode the Variant_effect column and other columns if needed 
def explode_columns_aligned(df, columns_to_explode):
    # Split each of the specified columns by commas
    for column in columns_to_explode:
        df[column] = df[column].apply(lambda x: x.split(',') if isinstance(x, str) else [x])

    # Apply exploded in a way that ensures all columns explode together
    df_exploded = df.copy()
    
    # Explode the rows based on the first column
    df_exploded = df_exploded.explode(columns_to_explode[0], ignore_index=True)
    
    # Now explode the other columns one by one, using the same index
    for column in columns_to_explode[1:]:
        df_exploded[column] = df_exploded[column].explode(ignore_index=True)
    
    return df_exploded


# Helper function which extracts the Transcript_ID from variant effect´s
def extract_gene_ids_from_variant_effect(variant_effect):
    # Extract gene ids (make sure this does not modify variant_effect directly). here you can add different patterns for plants 
    gene_id_pattern = r'(AT\dG\d+\.\d+|Os\d+t\d+(-\d+)?)'
    matches = re.findall(gene_id_pattern, variant_effect)
    
    return ', '.join([match[0] for match in matches if match[0]])  # Only join the gene IDs, ignoring empty strings



# Helper function which removes the variant effect´s which are not found inside the cds 
def remove_non_coding_terms(df):
    # List of non-coding region related SO terms to remove
    non_coding_terms = [
        "5_prime_UTR_variant", 
        "3_prime_UTR_variant", 
        "non_coding_transcript_exon_variant",
        "intron_variant", 
        "non_coding_transcript_variant", 
        "regulatory_region_variant", 
        "intergenic_variant", 
        "upstream_gene_variant", 
        "downstream_gene_variant",
        "mature_miRNA_variant"
        "TF_binding_site_variant",
        "TFBS_amplification"
        "NMD_transcript_variant"
        "TFBS_ablation"
        "regulatory_region_ablation"
        "regulatory_region_amplification"
    ]
    
    # Create a regular expression pattern to match any of the non-coding terms
    non_coding_pattern = '|'.join(non_coding_terms)
    
    # Remove the non-coding terms from the 'Variant_effect' column
    df['Variant_effect'] = df['Variant_effect'].apply(lambda x: ', '.join([term for term in x.split(', ') if not any(non_coding_term in term for non_coding_term in non_coding_terms)]) if isinstance(x, str) else x)
    
    # drop rows where the 'Variant_effect' column is empty after the removal
    df = df[df['Variant_effect'].str.strip() != '']

    return df


# Helper function to retrieve cds and pep sequences for multiple gene IDs in one row
def get_sequences_for_gene_ids(gene_ids, cds_dict, pep_dict):

    # Split multiple gene IDs if they exist
    gene_ids = gene_ids.split(",") 

    # CDS sequences, "NA" if not found
    gene_seqs = [cds_dict.get(gene_id.strip(), "NA") for gene_id in gene_ids]

    # Protein sequences, "NA" if not found
    pep_seqs = [pep_dict.get(gene_id.strip(), "NA") for gene_id in gene_ids]  
    
    return pd.Series([', '.join(gene_seqs), ', '.join(pep_seqs)])


# Helper function to parse protein (PEP) sequences
def parse_pep(fasta_file):
    return {record.id.split()[0]: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}


# Helper function to parse CDS sequences and extract start and end positions from the FASTA header. also the Start and end position in the chromosome if needed
# Compile can be updated so this function works on other plants
def parse_fasta(fasta_file):
    seq_dict = {}
    start_end_dict = {}

    # Precompile regex for speed
    regex = re.compile(r'chromosome:TAIR10:\d+:(\d+):(\d+):([+-1])')  

    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            gene_id = record.id.split()[0]
            seq_dict[gene_id] = str(record.seq)

            match = regex.search(record.description)
            if match:
                start_position, end_position, _ = match.groups()
                start_end_dict[gene_id] = (int(start_position), int(end_position))

    return seq_dict, start_end_dict


# Helper function which calculates the position of mutations if the mutation occurs inside the cds 
def adjust_mutation_position(mutation_start, CDSs, introns, strand):
    if len(CDSs) > 1:
        # Identify the start of the first CDS and end of the last CDS (same for both strands)
        first_CDS_start = CDSs.iloc[0].Start  # The first CDS (sorted by Start)

        if strand == "+":
            # Calculate the mutation position relative to the CDS (same as before for the positive strand)
            adjusted_position = mutation_start - first_CDS_start

            # Process the mutation for the positive strand
            for CDS in CDSs.itertuples():
                if CDS.Start <= mutation_start <= CDS.End:
                    # Get all introns before this CDS
                    intron_length_before_CDS = 0
                    for intron in introns.itertuples():
                        if intron.Start < CDS.Start:
                            intron_length_before_CDS += (intron.End - intron.Start + 1)
                        else:
                            break  # Stop once we encounter an intron after the CDS

                    # Adjust the mutation position by subtracting the length of the introns before this CDS
                    adjusted_position = (adjusted_position - intron_length_before_CDS) + 1
                    return int(adjusted_position)  # Ensure integer return value

            return None  # Return None if no CDS contains the mutation

        elif strand == "-":
            # Calculate the mutation position relative to the CDS (reverse calculation for the negative strand)
            last_CDS_end = CDSs.iloc[-1].End  # The last CDS (sorted by End)
            adjusted_position = last_CDS_end - mutation_start  # Reverse orientation calculation

            for CDS in CDSs.itertuples():  # Iterate through CDSs in order
                if CDS.Start <= mutation_start <= CDS.End:
                    # We found the CDS with the mutation
                    # Get all introns after this CDS (introns after the mutation for negative strand)
                    intron_length_after_CDS = 0
                    for intron in introns.itertuples():  # Iterate through introns in order
                        if intron.Start > CDS.End:  # Intron after the CDS (Start after CDS.End)
                            intron_length_after_CDS += (intron.End - intron.Start + 1)
                        else:
                            continue  # Continue iterating, do not break early. The first intron after the CDS is what matters.

                    # Adjust the mutation position by subtracting the length of the introns after this CDS
                    adjusted_position = (adjusted_position - intron_length_after_CDS) + 1
                    return int(adjusted_position)  # Ensure integer return value

            return None  # Return None if no CDS contains the mutation

        return None  # If strand is not '+' or '-', return None


# Helper function which gets the amino acid at the mutation position 
def get_amino_acid_at_position(pep_dict, transcript_id, adjusted_mut_position):
    # Divide by 3 and round up at the start to get the amino acid position
    adjusted_aa_position = math.ceil(adjusted_mut_position / 3)
    
    # Get the peptide sequence for the given transcript_id
    pep_sequence = pep_dict.get(transcript_id, "")
    
    # Check if the peptide sequence exists and if the adjusted position is valid
    if pep_sequence and 1 <= adjusted_aa_position <= len(pep_sequence):
        return pep_sequence[adjusted_aa_position - 1]  # Adjust for 0-based indexing
    return "NA"


# Helper function which gets the base at the mutation position and gives the complementary base for the minus strand 
def get_base_at_position(cds_dict, transcript_id, adjusted_mut_position, strand):
    cds_sequence = cds_dict.get(transcript_id, "")
    
    if cds_sequence and 1 <= adjusted_mut_position <= len(cds_sequence):
        base = cds_sequence[adjusted_mut_position - 1]  # Adjust for 0-based indexing
        
        # If the strand is negative, return the complementary base
        if strand == '-':
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            base = complement.get(base, base)  # Default to the base itself if it's not A, T, C, or G
            
        return base
    
    return "NA"  # Return "NA" if no valid base is found



# Helper function to extract transcript_id from the Suppl column in the intron and cds file
def extract_transcript_id(attribute_string):
    match = re.search(r'transcript_id "([^"]+)"', attribute_string)
    if match:
        return match.group(1)
    return None


# Helper function which assigns the Impact based on the variant effect of the mutation 
def assign_impact(variant_effect_cleaned):
    # Define the impact_mapping inside the helper function
    impact_mapping = {
        'transcript_ablation': 'HIGH',
        'splice_acceptor_variant': 'HIGH',
        'splice_donor_variant': 'HIGH',
        'stop_gained': 'HIGH',
        'frameshift_variant': 'HIGH',
        'stop_lost': 'HIGH',
        'start_lost': 'HIGH',
        'transcript_amplification': 'HIGH',
        'feature_elongation': 'HIGH',
        'feature_truncation': 'HIGH',
        'inframe_insertion': 'MODERATE',
        'inframe_deletion': 'MODERATE',
        'missense_variant': 'MODERATE',
        'protein_altering_variant': 'MODERATE',
        'splice_donor_5th_base_variant': 'LOW',
        'splice_region_variant': 'LOW',
        'splice_donor_region_variant': 'LOW',
        'splice_polypyrimidine_tract_variant': 'LOW',
        'incomplete_terminal_codon_variant': 'LOW',
        'start_retained_variant': 'LOW',
        'stop_retained_variant': 'LOW',
        'synonymous_variant': 'LOW',
        'coding_sequence_variant': 'MODIFIER',
        'mature_miRNA_variant': 'MODIFIER',
        '5_prime_utr_variant': 'MODIFIER',
        '3_prime_utr_variant': 'MODIFIER',
        'non_coding_transcript_exon_variant': 'MODIFIER',
        'intron_variant': 'MODIFIER',
        'nmd_transcript_variant': 'MODIFIER',
        'non_coding_transcript_variant': 'MODIFIER',
        'coding_transcript_variant': 'MODIFIER',
        'upstream_gene_variant': 'MODIFIER',
        'downstream_gene_variant': 'MODIFIER',
        'tfbs_ablation': 'MODIFIER',
        'tfbs_amplification': 'MODIFIER',
        'tf_binding_site_variant': 'MODIFIER',
        'regulatory_region_ablation': 'MODIFIER',
        'regulatory_region_amplification': 'MODIFIER',
        'regulatory_region_variant': 'MODIFIER',
        'intergenic_variant': 'MODIFIER',
        'sequence_variant': 'MODIFIER'
    }

    if variant_effect_cleaned:
        return ', '.join([impact_mapping.get(effect.strip(), 'UNKNOWN') for effect in variant_effect_cleaned.split(',')])
    return ''



start_time = time.time()
column_names = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']
df = pd.read_csv(gvf_file, sep='\t', index_col=False, names=column_names, dtype={'Start': int, 'End': int}, comment='#')
df = df.drop(columns=['Strand'])

# Define patterns for each field we want to extract
patterns = {
    'ID': r'ID=([^;]+)',
    'Reference_seq': r'Reference_seq=([^;]+)',
    'Variant_seq': r'Variant_seq=([^;]+)',
    'dbxref': r'Dbxref=([^;]+)',
    'reference_peptide': r'reference_peptide=([^;]+)',
    'variant_peptide': r'variant_peptide=([^;]+)',
    'sift_prediction': r'sift_prediction=([^;]+)'
}

# Define patterns for Variant effect first so it can be exploded 
variant_effect_extract = {'Variant_effect': r'Variant_effect=([^;]+)'}

# Step 1: Extract 'Variant_effect' from the 'Suppl' column and immediately remove it from 'Suppl'
df['Variant_effect'] = df['Suppl'].apply(lambda x: extract_id(x, variant_effect_extract['Variant_effect']))

# Remove 'Variant_effect' from 'Suppl' once it's extracted
df['Suppl'] = df['Suppl'].apply(lambda x: re.sub(r'Variant_effect=[^;]+;?', '', x).strip(';'))

print(f"Step 1 finished at: {time.time() - start_time:.2f} seconds")

# Step 2: Drop rows where 'Variant_effect' is NaN
df = df.dropna(subset=['Variant_effect'])

print(f"Step 2 finished at: {time.time() - start_time:.2f} seconds")

# Step 3: Apply the explode-function to each of the specified columns
columns_to_explode = ['Variant_effect']
df = explode_columns_aligned(df=df, columns_to_explode=columns_to_explode)

print(f"Step 3 finished at: {time.time() - start_time:.2f} seconds")

# Step 4: Remove non-coding terms from 'Variant_effect'
df = remove_non_coding_terms(df)

print(f"Step 4 finished at: {time.time() - start_time:.2f} seconds")

# Step 5: Add Transcript_ID derived from Variant_effect. here you can add different gene_id_pattern´s in the extract_gene_ids_from_variant_effect funktion for other plants
df['Transcript_ID'] = df['Variant_effect'].apply(extract_gene_ids_from_variant_effect)

print(f"Step 5 finished at: {time.time() - start_time:.2f} seconds")

# Step 6: Extract the other field from Suppl and create new columns
for field, pattern in patterns.items():
    df[field] = df['Suppl'].apply(lambda x: extract_id(x, pattern))

# Remove extracted fields from 'Suppl' and keep only the remaining parts
remove_pattern = r'(ID|Reference_seq|Variant_seq|Dbxref|reference_peptide|variant_peptide|sift_prediction)=[^;]+;?'
df['Suppl'] = df['Suppl'].apply(lambda x: re.sub(remove_pattern, '', x).strip(';'))

print(f"Step 6 finished at: {time.time() - start_time:.2f} seconds")

# Step 7: Deletes all variant peptides except the one with matching Transcript_ID and drop rows were variant_peptide is empty
df = df.dropna(subset=['variant_peptide'])

df = filter_variant_peptide_by_transcript(df)

print(f"Step 7 finished at: {time.time() - start_time:.2f} seconds")

# Step 8: Clean the Variant_effect column (remove transcript info and normalize)
df['Variant_effect_cleaned'] = df['Variant_effect'].apply(clean_variant_effect)

print(f"Step 8 finished at: {time.time() - start_time:.2f} seconds")

# Step 9: Assign the IMPACT value based on the cleaned Variant_effect
df['Impact'] = df['Variant_effect_cleaned'].apply(assign_impact)

print(f"Step 9 finished at: {time.time() - start_time:.2f} seconds")

# Step 10: Parse CDS and PEP files
cds_dict, cds_start_end_dict = parse_fasta(cds_file)  # Parse CDS and extract start/end
pep_dict = parse_pep(pep_file)  # Parse PEP files

print(f"Step 10 finished at: {time.time() - start_time:.2f} seconds")

# Step 11: Remove rows where all columns except 'ID' are the same
df = df.loc[df.drop('ID', axis=1).duplicated(keep='first') == False]

print(f"Step 11 finished at: {time.time() - start_time:.2f} seconds")

# Step 12: Get sequences for each row (CDS and PEP). Can be excluded if not needed to reduce the size of the created file but need to be removed from column_order 
df[['Gen_Sequence', 'Amino_Acid_Sequence']] = df['Transcript_ID'].apply(get_sequences_for_gene_ids, cds_dict=cds_dict, pep_dict=pep_dict)

print(f"Step 12 finished at: {time.time() - start_time:.2f} seconds")

# Step 13: Strip any leading/trailing whitespace from the 'Transcript_ID' column
df['Transcript_ID'] = df['Transcript_ID'].str.strip()

print(f"Step 13 finished at: {time.time() - start_time:.2f} seconds")


# Step 14: Read the cds intron gtf file
CDS_intron_df = add_introns(gtf_file)
print(f"Step 14 finished at: {time.time() - start_time:.2f} seconds")

# Step 15: Extract transcript_id from the Suppl column in the intron and cds file
CDS_intron_df['transcript_id'] = CDS_intron_df['Suppl'].apply(extract_transcript_id).str.upper().str.strip()

# Create a new DataFrame with 'transcript_id' and 'Strand'
transcript_strand = CDS_intron_df[['transcript_id', 'Strand']].drop_duplicates()

print(f"Step 15 finished at: {time.time() - start_time:.2f} seconds")

# Step 16: Initialize the list for adjusted mutation positions
adjusted_mut_positions = []

print(f"Step 16 finished at: {time.time() - start_time:.2f} seconds")


# Step 17: Process mutations in the GVF file
# Precompute CDS/intron data as a dictionary for O(1) lookups
CDS_intron_dict = {
    (tid.upper(), sid): grp for (tid, sid), grp in CDS_intron_df.groupby(['transcript_id', 'Seqid'])
}

# Precompute transcript strand as a dictionary for O(1) lookups
transcript_strand_dict = dict(zip(transcript_strand['transcript_id'].str.upper(), transcript_strand['Strand']))

rows_to_drop = []
adjusted_mut_positions = []
strand_updates = {}

# Iterate efficiently using itertuples()
for row in df.itertuples(index=True):
    idx = row.Index
    seqid = row.Seqid
    mutation_start = row.Start
    transcript_id = row.Transcript_ID.upper()

    # Fast lookup for CDSs/introns
    CDS_introns = CDS_intron_dict.get((transcript_id, seqid), None)
    if CDS_introns is None:
        rows_to_drop.append(idx)
        continue

    # Extract and sort CDSs and introns
    CDSs = CDS_introns[CDS_introns['Type'] == 'CDS'].sort_values(by='Start')
    introns = CDS_introns[CDS_introns['Type'] == 'intron'].sort_values(by='Start')

    # Fast lookup for strand
    strand = transcript_strand_dict.get(transcript_id, None)
    if strand:
        strand_updates[idx] = strand

    # Adjust mutation position
    adjusted_mut_position = adjust_mutation_position(mutation_start, CDSs, introns, strand)
    adjusted_mut_positions.append(adjusted_mut_position)

# Bulk update strand column
df.loc[strand_updates.keys(), 'Strand'] = list(strand_updates.values())

# Drop rows in bulk
df.drop(rows_to_drop, inplace=True)
df['Adjusted_Mutation_Position'] = adjusted_mut_positions

print(f"Step 17 finished at: {time.time() - start_time:.2f} seconds")


# Step 18: Add the adjusted mutation positions as a new column in the DataFrame and drop all rows were the Mutation is not inside an CDS 
df = df.dropna(subset=["Adjusted_Mutation_Position"])

df["Adjusted_Mutation_Position"] = df["Adjusted_Mutation_Position"].astype(int)

print(f"Step 18 finished at: {time.time() - start_time:.2f} seconds")

# Step 19: Add Base_at_Adjusted_Position as a new column in the DataFrame
df['Base_at_Adjusted_Position'] = df.apply(
    lambda row: get_base_at_position(cds_dict, row['Transcript_ID'], row['Adjusted_Mutation_Position'], row['Strand']), axis=1)

print(f"Step 19 finished at: {time.time() - start_time:.2f} seconds")



# Step 20: Add AA_at_Adjusted_Position as a new column in the DataFrame
df['AA_at_Adjusted_Position'] = df.apply(
    lambda row: get_amino_acid_at_position(pep_dict, row['Transcript_ID'], row['Adjusted_Mutation_Position']), axis=1)

print(f"Step 20 finished at: {time.time() - start_time:.2f} seconds")

# Step 21: Drop Columns which are not needed
df = df.drop(columns=['Source', 'Score', 'Frame', 'Suppl', 'sift_prediction', 'Variant_effect_cleaned'])

print(f"Step 21 finished at: {time.time() - start_time:.2f} seconds")

# Step 22: Order Columns
colum_order = ['Seqid', 'ID', 'Type', 'Start', 'End', 'Adjusted_Mutation_Position', 'Reference_seq', 'Base_at_Adjusted_Position', 'Variant_seq', 'reference_peptide', 'AA_at_Adjusted_Position', 'Strand', 'variant_peptide', 'Transcript_ID', 'Variant_effect', 'Impact', 'dbxref', 'Gen_Sequence', 'Amino_Acid_Sequence']
df = df[colum_order]

print(f"Step 22 finished at: {time.time() - start_time:.2f} seconds")

# Step 23: Sort by Seqid and then by Start Value
df = df.sort_values(by=['Seqid', 'Start'], ascending=[True, True])

print(f"Step 23 finished at: {time.time() - start_time:.2f} seconds")

# Step 24: Save the Final Dataframe
base_filename = os.path.splitext(gvf_file)[0]
output_file = base_filename + '_main.gvf'
df.to_csv(output_file, sep='\t', header=True, index=False)

print(f"Step 24 finished at: {time.time() - start_time:.2f} seconds")


