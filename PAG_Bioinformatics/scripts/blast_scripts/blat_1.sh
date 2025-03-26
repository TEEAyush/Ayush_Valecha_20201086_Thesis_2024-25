#!/bin/bash
#SBATCH --job-name=blat
#SBATCH --partition=COMPUTE
#SBATCH --output=blat_1_out.txt
#SBATCH --error=blat_1_error.txt
#SBATCH --time=05:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4


module load blat/36

working_dir="/home/ayush.valecha/workingdir/pollen_allergen_genes/data"
GENOME="${working_dir}/sequence_files/Pan_040_2CCS.HiFiasm.primary_renamed.fasta"
CDS_SEQUENCES="${working_dir}/sequence_files/PAG_seqkit_exons_nseq.fasta"
OUTPUT="${working_dir}/sequence_files/blat_1_out.psl"

blat $GENOME $CDS_SEQUENCES $OUTPUT
