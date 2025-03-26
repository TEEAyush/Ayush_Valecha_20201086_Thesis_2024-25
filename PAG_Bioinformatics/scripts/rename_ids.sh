#!/bin/bash

work_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data" 
ID_FILE="${work_dir}/PAG_IDS.txt"
#FASTA_FILES=("${work_dir}/sequence_files/PAG_nseq.fasta" "${work_dir}/sequence_files/PAG_nseq_s.fasta" "${work_dir}/sequence_files/PAG_protein.fasta" "${work_dir}/sequence_files/PAG_seqkit_exons_nseq.fasta" )
FASTA_FILES=("${work_dir}/sequence_files/not-diff-expr_PAG_exons_nseq.fasta")
OUTPUT_DIR="${work_dir}/sequence_files/renamed_PAG"
MAP_FILE="${work_dir}/renaming_PAG.txt"

# Create output directory if it doesn't exist
#mkdir -p "$OUTPUT_DIR"

# Create the renaming map
#awk '{print $1, "PAG" NR}' "$ID_FILE" > "$MAP_FILE"

# Function to rename IDs in FASTA file
rename_fasta() {
	    local fasta_file=$1
	    local output_file=$2
	    cp "$fasta_file" "$output_file"
	    while read -r old_id new_id; do
		    sed -i -e "s/\(^>\)$old_id\([[:alnum:]_]*\)/\1$new_id\2/g" "$output_file"
	    done < "$MAP_FILE"
}

# Rename IDs in each FASTA file
for FASTA in "${FASTA_FILES[@]}"; do
	OUTPUT_FILE="$OUTPUT_DIR/$(basename "$FASTA")"
	rename_fasta "$FASTA" "$OUTPUT_FILE"
done

echo "Renaming completed. Renamed FASTA files are in the '$OUTPUT_DIR' directory."
#echo "Renaming map saved to '$MAP_FILE'."

