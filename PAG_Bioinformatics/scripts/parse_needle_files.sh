#!/bin/bash


work_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/sequence_files/renamed_PAG"
# Output directory for needle results
OUTPUT_DIR="$work_dir/needle_results"
# Output file for similarity scores
SIMILARITY_FILE="$work_dir/similarity_scores.txt"

# Initialize the similarity scores file with a header
echo -e "Seq1\tSeq2\tIdentity\tSimilarity\tGaps\tScore" > "$SIMILARITY_FILE"

# Parse each needle output file to extract metrics
for needle_output in "$OUTPUT_DIR"/*.needle; do
	   seq1=$(basename "$needle_output" | cut -d'_' -f1)
	   seq2=$(basename "$needle_output" | cut -d'_' -f2 | cut -d'.' -f1)

	   dentity=$(grep -m 1 "Identity:" "$needle_output" | awk '{print $2}' | tr -d '()%' | awk -F/ '{print $1/$2}')
	   similarity=$(grep -m 1 "Similarity:" "$needle_output" | awk '{print $2}' | tr -d '()%' | awk -F/ '{print $1/$2}')
	   gaps=$(grep -m 1 "Gaps:" "$needle_output" | awk '{print $2}' | tr -d '()%')
	   score=$(grep -m 1 "Score:" "$needle_output" | awk '{print $2}')

	   echo -e "${seq1}\t${seq2}\t${identity}\t${similarity}\t${gaps}\t${score}" >> "$SIMILARITY_FILE"
done

echo "Similarity scores calculated and saved to '$SIMILARITY_FILE'."

