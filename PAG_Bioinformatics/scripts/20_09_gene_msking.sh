#!/bin/bash

# Input files
work_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/sequence_files"
assembly_to_mask="$work_dir/hap1/herma_files/Pan_01N25.hap1_rt.fasta"
gene_sequences="$work_dir/PAG_full_gene_seq.fasta"

# Output file
masked_assembly="$work_dir/hap1/herma_files/Pan_01N25.hap1_rt_mask.fasta"

# Copy assembly B to the output file
cp "$assembly_to_mask" "$masked_assembly"

# Loop through each gene in PAG_full_gene_seq.fasta
while read -r header; do
	  if [[ $header == ">"* ]]; then
		gene_name="$header"
	  else
		gene_seq="$header"
	  # Create a string of Ns with the same length as the gene sequence
	        mask_seq=$(echo "$gene_seq" | sed 's/./N/g')

	        # Use sed to replace the gene sequence with Ns in the masked assembly file
	        sed -i "s/$gene_seq/$mask_seq/g" "$masked_assembly"
	  fi
done < <(grep -A 1 ">" "$gene_sequences")

