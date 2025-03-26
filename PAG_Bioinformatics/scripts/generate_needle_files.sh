#!/bin/bash

work_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/sequence_files/renamed_PAG"
# Input FASTA file
FASTA_FILE="${work_dir}/not-diff-expr_PAG_exons_nseq.fasta"
# Output directory for needle results
OUTPUT_DIR="${work_dir}/needle_results"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Extract sequence IDs and sequences from the FASTA file
awk '/^>/ {if (seq) print seq; print; seq=""} {if (!/^>/) seq=seq$0} END {print seq}' "$FASTA_FILE" > "${work_dir}/temp.fasta"

# Read sequence IDs into an array
readarray -t SEQ_IDS < <(grep '^>' "${work_dir}/temp.fasta" | tr -d '>')

#Run needle for each pair of sequences
for ((i=0; i<${#SEQ_IDS[@]}-1; i++)); do
	for ((j=i+1; j<${#SEQ_IDS[@]}; j++)); do
		 seq1="${SEQ_IDS[$i]}"
	         seq2="${SEQ_IDS[$j]}"
		
		 seq1_file="$OUTPUT_DIR/${seq1}.fasta"
        	 seq2_file="$OUTPUT_DIR/${seq2}.fasta"

        	 # Extract sequences for the current pair
        	 awk -v id="$seq1" '/^>/ {print; getline; if ($0 ~ id) print $0}' "${work_dir}/temp.fasta" > "$seq1_file"
        	 awk -v id="$seq2" '/^>/ {print; getline; if ($0 ~ id) print $0}' "${work_dir}/temp.fasta" > "$seq2_file"
        	
		 
		 needle_output="$OUTPUT_DIR/${seq1}_${seq2}.needle"


					             
		 # Run needle for the current pair of sequences
		 needle -asequence "$seq1_file" -bsequence "$seq2_file" -gapopen 10 -gapextend 0.5 -outfile "$needle_output" \
			-aformat3 pair 

		rm "$seq1_file" "$seq2_file" 
	done
done

# Clean up temporary files
rm "${work_dir}/temp.fasta"

echo "Needle files generated in '$OUTPUT_DIR'."
