#!/bin/bash

work_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/sequence_files"
# Input GFF file
#input_files="Pan_040.hap1_rt_no_head.gff" "Pan_040.hap2_rt_no_head.gff" "Pan_01N25.hap1_rt_no_head.gff" "Pan_01N25.hap2_rt_no_head.gff"

input_files="${work_dir}/blat_PAG_mask_renamed_nh.gff"

for GFF_FILE in $input_files; do
	base_name=$(basename "$GFF_FILE" .gff)

	OUTPUT_FILE="${work_dir}/filtered_${base_name}.gff"

#	TEMP_FILE="${work_dir}/temp_filtered_${base_name}.gff"

# Step 1: Filter by chromosome name (chr_M_i)
#awk '$1 ~ /^chr_M_[0-9]+$/ {print $0}' "$GFF_FILE" > "$TEMP_FILE"

	# Step 2: Group exons by transcript ID and calculate exon count, average exon length, and gene length
	awk '
	BEGIN {
	    FS = OFS = "\t"
	}
	$3 == "exon" {
	    # Extract transcript ID
	    match($9, /Parent=([^;]+)/, arr)
	    transcript_id = arr[1]
	    chrom = $1
	    start = $4
	    end = $5
	    strand = $7
	
	    exon_lengths[transcript_id] += end - start + 1
	    exon_counts[transcript_id]++
	    genes[transcript_id, "start"] = (genes[transcript_id, "start"] == "" || genes[transcript_id, "start"] > start) ? start : genes[transcript_id, "start"]
	    genes[transcript_id, "end"] = (genes[transcript_id, "end"] == "" || genes[transcript_id, "end"] < end) ? end : genes[transcript_id, "end"]
	}
	END {
	     for (transcript_id in exon_counts) {
		avg_exon_length = exon_lengths[transcript_id] / exon_counts[transcript_id]
		gene_length = genes[transcript_id, "end"] - genes[transcript_id, "start"] + 1
		if (exon_counts[transcript_id] == 2 || exon_counts[transcript_id] == 3) {
		   if (avg_exon_length >= 200 && avg_exon_length <= 400) {
	              if (gene_length >= 900 && gene_length <= 1200) {
			 print transcript_id
			}
		   }
		}
	     }
	}' "$GFF_FILE" > transcript_ids.txt

	# Step 3: Filter original GFF file to include only the desired transcripts
	awk 'FNR==NR {transcripts[$1]; next} 
	{
	     match($9, /Parent=([^;]+)/, arr)
	     transcript_id = arr[1]
	     if (transcript_id in transcripts) {
		print $0
	     }
	}' transcript_ids.txt "$GFF_FILE" > "$OUTPUT_FILE"

	# Clean up temporary files
	rm transcript_ids.txt

	echo "Filtering completed. Filtered GFF file saved to '$OUTPUT_FILE'."

done
