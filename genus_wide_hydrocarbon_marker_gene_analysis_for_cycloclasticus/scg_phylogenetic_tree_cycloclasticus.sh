#1. Using GtoTree workflow: Annotate single copy core genes, make a multiple sequence alignment,
#create a maximum likelihood tree with IQTree, and model finder.
##We will use the HMM set for Gammaproteobacteria provided by GtoTree to annotate the scgs
##Written as a slurm array job for a cluster on Compute Canada (Digital Research Alliance of Canada).

#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --account=def-nagissa
#SBATCH --job-name=gtotree

set -euo pipefail
set -x

module load apptainer

#Required inputs
SIF=/project/6105702/general_scripts_and_sifs/gtotree.sif

#I realize these variables are annoying- but its set up to make sure the paths are
#bound properly to apptainer

#File legend (realpaths to all fasta files in one file)

file_legend_path=/project/6105702/CMO_2024_Metagenome/CYCLOCLASTICUS_PHYLOGENOMICS
file_legend_file_name=fasta_paths.txt

#Add path to directory where fasta files are stored
fasta_files_path=/project/6105702/CMO_2024_Metagenome/CYCLOCLASTICUS_PHYLOGENOMICS/fasta

#Make output directory
outdir=/project/6105702/CMO_2024_Metagenome/CYCLOCLASTICUS_PHYLOGENOMICS/scgs
mkdir -p ${outdir}

#Run gtotree
apptainer exec \
  --bind ${file_legend_path} \
  --bind ${fasta_files_path} \
  --bind ${outdir} \
  ${SIF} \
  GToTree \
        -f ${file_legend_path}/${file_legend_file_name} \
        -H Gammaproteobacteria \
        -j 4 \
        -o ${outdir} \
        -T IQ-TREE \
        -m MFP \
        -B 1000 \
        -pre iqtree_cyclo_gamma_scg \
        -n $SLURM_CPUS_PER_TASK \
        -F #THIS WILL OVERRIDE THE OUTPUT DIRECTORY IF IT ALREADY EXISTS
