#!/bin/bash

#file paths as variables
working_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data/coordinate_files"

rm ${working_dir}/Ha_Pan_040_marked.gff ${working_dir}/Hb_Pan_040_marked.gff

touch ${working_dir}/Ha_Pan_040_marked.gff ${working_dir}/Hb_Pan_040_marked.gff

for SIgrp in Ha Hb; do
	rm ${working_dir}/${SIgrp}_Pan_040_marked.gff

	touch ${working_dir}/${SIgrp}_Pan_040_marked.gff

	grep -F -f ${working_dir}/${SIgrp}_Pan_040_PAG_transcript.gff ${working_dir}/${SIgrp}_Pan_040_StringTie_transcript.gff | awk '{print $0"\t"1}'q >>  ${working_dir}/${SIgrp}_Pan_040_marked.gff

	grep -v -F -f ${working_dir}/${SIgrp}_Pan_040_PAG_transcript.gff ${working_dir}/${SIgrp}_Pan_040_StringTie_transcript.gff | awk '{print $0"\t"0}' >>  ${working_dir}/${SIgrp}_Pan_040_marked.gff

	grep -F -f ${working_dir}/${SIgrp}_Pan_040_PAG_exon.gff ${working_dir}/${SIgrp}_Pan_040_StringTie_exon.gff | awk '{print $0"\t"1}' >>  ${working_dir}/${SIgrp}_Pan_040_marked.gff

	grep -v -F -f ${working_dir}/${SIgrp}_Pan_040_PAG_exon.gff ${working_dir}/${SIgrp}_Pan_040_StringTie_exon.gff | awk '{print $0"\t"0}' >>  ${working_dir}/${SIgrp}_Pan_040_marked.gff

	cat ${working_dir}/${SIgrp}_Pan_040_maker_exon.gff
