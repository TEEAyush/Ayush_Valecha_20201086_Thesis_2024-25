#!/bin/bash

work_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/sequence_files"
# Input GFF file
input_gff="${work_dir}/filtered_blat_PAG_mask_renamed_nh.gff"
# Output merged GFF file
base_name=$(basename "$input_gff" .gff)
output_gff="${work_dir}/${base_name}_uniq.gff"

# Convert GFF to BED-like format
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $9}' "$input_gff" > "${work_dir}/temp.bed"

# sort
sort -k1,1 -k2,2n "${work_dir}/temp.bed" > "${work_dir}/temp1.bed"

# Merge overlapping entries while preserving the attributes
bedtools merge -i "${work_dir}/temp1.bed" -c 4 -o distinct > "${work_dir}/merged.bed"

# Convert merged BED-like format back to GFF
awk 'BEGIN{OFS="\t"} {print $1, ".", "exon", $2+1, $3, ".", ".", ".", $4}' "${work_dir}/merged.bed" > "$output_gff"

# Clean up temporary files
rm "${work_dir}/temp.bed" "${work_dir}/temp1.bed" "${work_dir}/merged.bed"

echo "Merging completed. Merged GFF file saved to '$output_gff'."
