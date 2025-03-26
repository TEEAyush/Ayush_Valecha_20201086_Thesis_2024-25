#!/bin/bash

work_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/sequence_files/hap1"
input="${work_dir}/Pan_040_hap1_all_cds.fasta"
output="${work_dir}/Pan_040_hap1_all_cds1.fasta"





awk '{
    if (substr($0, 1, 1) == ">") { 
	     print $0;  # Print the header as is
    } else { 
    if (length($0) > 3 && substr($0, 2, 3) == "ATG") {
	     print substr($0, 2);  # Remove the first base if it precedes "ATG"
    } else {
             print $0;  # Otherwise, print the sequence as is
    }
   }
}' "$input" > "$output"

