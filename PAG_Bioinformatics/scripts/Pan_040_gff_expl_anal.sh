#!/bin/bash

#file paths as variables
gff_file="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/Pan_040_primary.maker.output1/Pan_040_primary.all.gff"
output_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/results_Pan_040_gff_expl_anal"

# Create the type column list
cat ${gff_file}  | cut -f3  | sort | uniq | awk '$0 !~ "\#"'  > ${output_dir}/Pan_040_element_type.txt   # don't want the lines with comments

# iterate over the list, create a file with lengths for each type 
for type in $(cat ${output_dir}/Pan_040_element_type.txt) ; do 
	cat ${gff_file} | \
	awk -v type=$type '{if($3==type) print $5-$4"\t"$9}' \
	> ${output_dir}/Pan_040_element_${type}_lengths.txt;
        
done




