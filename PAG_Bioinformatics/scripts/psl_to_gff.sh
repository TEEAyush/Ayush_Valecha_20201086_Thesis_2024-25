#!/bin/bash

work_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/sequence_files/hap1"
PSL_FILE="${work_dir}/blat_Pan_040_filtered_20_09.psl"


base_name=$(basename "$PSL_FILE" .psl)
GFF_FILE="${work_dir}/${base_name}.gff"
	

# Initialize alignment number
aln_no=1


{
# Skip header lines
grep -v -e '^psLayout' -e '^---' "$PSL_FILE" | awk -v aln_no="$aln_no" '
BEGIN {
       print "##gff-version 3"
      }
      {
       chrom = $14
       strand = $9
       name = $10
       block_sizes = $19
       block_starts = $21
       start = $16
       end = $17
       score = "."  # Placeholder for score, GFF field 6
       # transcript_id = "transcript_aln" aln_no
       feature = "exon"
       source = "BLAT"
       split(block_sizes, sizes, ",")
       split(block_starts, starts, ",")
       # Print the parent transcript feature
       # print chrom "\t" source "\ttranscript\t" start "\t" end "\t" score "\t" strand "\t" "." "\tID=" transcript_id ";Name=" name      
       for (i = 1; i <= length(sizes); i++) {
		                if (sizes[i] > 0) {  # Exclude blocks of 0 length
				block_start = starts[i]
				block_end = block_start + sizes[i]
				query_name = name "_aln_" aln_no
				print chrom "\t" source "\t" feature "\t" block_start "\t" block_end "\t" score "\t" strand "\t" "." "\tParent=" query_name
			}
		}
	aln_no++
	}' 
	} > "$GFF_FILE"

echo "GFF file created: $GFF_FILE"

