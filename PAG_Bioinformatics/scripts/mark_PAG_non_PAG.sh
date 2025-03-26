#!/bin/bash

#file paths as variables
output_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/results_Pan_040_gff_expl_anal"
PAG_IDS="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/PAG_IDS.txt"
#All_gene="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/Pan_040_gene_IDS.txt"
#Type="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/results_Pan_040_gff_expl_anal/Pan_040_element_type.txt"


# iterate over the type list, create a file with lengths for each type
for type in CDS exon five_prime_UTR gene mRNA three_prime_UTR ; do
		        cat ${output_dir}/Pan_040_element_${type}_lengths.txt | sed -E 's/ID=(.*gene-[0-9.]+)[^;]*;.*[^[:space:]]+/\1/' > ${output_dir}/Pan_040_element_${type}_lengths1.txt
			rm ${output_dir}/Pan_040_element_${type}_lengths.txt
		        mv ${output_dir}/Pan_040_element_${type}_lengths1.txt ${output_dir}/Pan_040_element_${type}_lengths.txt
			grep -F -f ${PAG_IDS} ${output_dir}/Pan_040_element_${type}_lengths.txt | sed "s/$/\t1/" > ${output_dir}/Pan_040_element_${type}_lengths_marked.txt
			grep -v -F -f ${PAG_IDS} ${output_dir}/Pan_040_element_${type}_lengths.txt | sed "s/$/\t0/" >> ${output_dir}/Pan_040_element_${type}_lengths_marked.txt
			rm ${output_dir}/Pan_040_element_${type}_lengths.txt
			mv ${output_dir}/Pan_040_element_${type}_lengths_marked.txt ${output_dir}/Pan_040_element_${type}_lengths.txt
done


