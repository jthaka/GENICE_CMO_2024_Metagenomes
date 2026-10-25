#Goal: organize hits for hydrocarbon degrading marker genes, extract relevant FNA and AA seqs,
#construct multiple sequence alignments, and filter alignments for quality control.
    ##Provides alignments for downstream analyses. There are numerous marker genes, so I wrote this script to go from an hmm output
    ##for the CANTHYD database to multiple sequence alignments for those genes across your mags/genomes of interest.

#Required inputs: 
    #a. List of genome names (needed to point tools to appropriate fasta or tsv files): all_cyclo_genomes.txt
    #b. Gene hits output from hmmer (in my case this was for the CANTHYD database)
    #c. List of gene calls (in my case I was working with the .gff3 output or .faa output from bakta
      #depending on if I was analysing nucleotide or amino acid sequences).
    #d. Python script: merge_hmmsearch_and_gene_annotation.py

#MULTIPLE BASH SCRIPTS/LOOPS WITH A COMMON SET OF BASH VARIABLES

#Set common variables
genome_list=/path/to/genome_list.txt  #This will set the ${genome} variable
genome_fasta_dir=path/to/genomes  #Genome fasta files referred to as ${genome_dir}/${genome}.fa in step 5b.
genome_fasta_file=${genome_dir}/${genome}.fa
hmm_tblout_dir=/path/to/hmm_hits
hmm_tblout_file=hmm_hits.tblout
gene_calls_dir=/path/to/gene/calls/output/per/${genome}  #NOTE: should have a separate dir for each genome
gene_calls_file=annotation.gff3
amino_acid_calls_dir=/path/to/gene/calls/output/per/${genome}
amino_acid_calls_file=annotation.fa
merge_python_script=path/to/merge_hmmsearch_and_gene_annotation.py
merge_python_script_requirements=path/to/merge_requirements.txt

#Make the primary output directory
outdir=path/to/primary/output/dir
mkdir -p ${outdir}

#1. Merge gff files with hmm outputs-----------------------------------------------------------------------------------------------------------------------------------------------
  ##(written as a slurm array job for the Nibi Compute Canada Cluster)

genome=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${genome_list})
chmod +x ${merge_python_script}
module load python/3.14
srun --ntasks $SLURM_NNODES --tasks-per-node=1 bash << EOF
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install -r ${merge_python_script_requirements}
EOF
# activate only on main node
source $SLURM_TMPDIR/env/bin/activate;
python ${script} \
        --gff3 ${gene_calls_dir}/${gene_calls_file} \
        --tblout ${hmm_tblout_dir}/${hmm_tblout_file} \
        --out ${outdir}/${genome}_hmm_w_gff.tsv \
        --join inner

#2a. Remove weird symbols from gene names if necessary such as spaces, '%', and '/', in gene names --------------------------------------------------------------------------------
  ##Helpful for organizing files by gene downstream. Below are examples I of characters I was removing.
  
for genome in $(cat ${genome_list}); do
hmm_w_gff=${outdir}/${genome}_hmm_w_gff.tsv
gawk -i inplace 'BEGIN{FS=OFS="\t"} {
gsub(/[[:blank:]]/, "_", $7)} 1' ${hmm_w_gff};
done
for genome in $(cat ${genome_list}); do
hmm_w_gff=${outdir}/${genome}_hmm_w_gff.tsv
gawk -i inplace 'BEGIN{FS=OFS="\t"} {
gsub(/'%'/, "_", $7)} 1' ${hmm_w_gff};
done
for genome in $(cat ${genome_list}); do
hmm_w_gff=${outdir}/${genome}_hmm_w_gff.tsv
gawk -i inplace 'BEGIN{FS=OFS="\t"} {
gsub("/", "_", $7)} 1' ${hmm_w_gff};
done

#2b. Specific for the CANTHYD output: I would like to create a field that combines the simple CANTHYD query name with gene description.--------------------------------------------
	##The gene description includes extra information of things like the particular subunit that the gene is encoding.
	##Therefore I want to have the option later on to organize by the detailed description of the hit.
	##I want to concat with the simple gene query, as some of the detailed descriptions can be identical between two enzyme systems
	##(e.g., MAH_alpha and TmoA could have the same gene product description, but I want to keep these separate to stay true to the HMMs annotation).

for genome in $(cat ${genome_list}); do
hmm_w_gff=${outdir}/${genome}_hmm_w_gff.tsv
gawk -i inplace 'BEGIN{FS=OFS="\t"} NR==1{$(NF+1)="enzyme_and_description"; print; next} {$(NF+1)="(" $3 ")" "_" $7; print}' ${hmm_w_gff};
done


#3. Concatenate files to get list of unique gene names-----------------------------------------------------------------------------------------------------------------------------

awk 'NR==1 || FNR>1' \
    ${outdir}/*_hmm_w_gff.tsv \
> ${outdir}/all_hmm_w_gff.tsv

##Make unique variable to keep referring to this unique list of HC-degrading marker genes
##Change $(NF) to whatever field you are organizing by.

awk 'BEGIN{FS=OFS="\t"} {print $(NF)}' ${outdir}/all_hmm_w_gff.tsv | cat | sort | uniq \
> ${outdir}/unique_hits.tsv

##Make unique variable to keep referring to this unique list of HC-degrading marker genes
unique=${outdir}/unique_hits.tsv


#4. Organize the genes into 1 type per file, but keep files separated by MAG/genome.-----------------------------------------------------------------------------------------------
##Note: column 3 referenced in the awk script below hosts the gene annotations that I am organizing by.

for genome in $(cat ${genome_list}); do
mkdir -p ${outdir}/${genome}
for gene in $(cat ${unique}); do
awk -F'\t' -v OFS='\t' \
    -v col=3 \
    -v val=${gene} \
    'NR==1 { print; next } $col==val { print }' \
${outdir}/${genome}_hmm_w_gff.tsv;
${outdir}/${genome}/${gene}_regions.tsv;
done;
done

#5a. Pull each amino acid sequence from the bakta .faa file using the bakta gene ID.-----------------------------------------------------------------------------------------------
    ##We are using a simple grep search to find the header and print the line immediately below. 
    ##We need to cycle through all unique bakta IDs in a given region file and pull out the amino acid sequences from bakta .faa output:
    ##(written as a slurm job for the Nibi Compute Canada Cluster)

for genome in $(cat ${genome_list}); do
for gene in $(cat ${unique}); do
for id in $(cat ${outdir}/${genome}/${gene}_regions.tsv | cut -f 1 | sort | uniq); do
grep -wA 1 "${id}" ${amino_acid_calls_dir}/${amino_acid_calls_file} \
>> ${outdir}/${genome}/${gene}_${genome}.faa;
done;
done;
done

#5b. Pull each genome sequence from the bakta .fna file using genome regions and the samtools faidx tool---------------------------------------------------------------------------
  ##(written as a slurm job for the Nibi Compute Canada Cluster)
  ##NOTE: in the awk script within faidx, change $NF to whatever column has the 
  
  ##create a region column that follows the format necessary for Samtools faidx (chr:from-to)
  ##Column number legend: chromosome/contig label from gene call software output ($9), start position ($11), end position ($12).

for genome in $(cat ${genome_list}); do
gawk -i inplace -F'\t' -v OFS='\t' -v genome="$genome" 'NR==1{$(NF+1)="samtools_faidx_region"; print; next} { $(NF+1) = $9 ":" $11 "-" $12; print }' \
${outdir}/${genome}/${gene}_regions.tsv;
done

module load samtools
for genome in $(cat ${genome_list}); do
for gene in $(cat ${unique}); do
samtools faidx -o ${outdir}/${genome}/${gene}_${genome}.fna \
-r <(awk -F'\t' 'NR>1{print $NF}' ${outdir}/${genome}/${gene}_regions.tsv) \
${genome_dir}/${genome}.fa;
done;
done


#6. In .fna and .faa files we just created, add genome name to beginning of the header.--------------------------------------------------------------------------------------------
for genome in $(cat ${genome_list}); do
for gene in $(cat ${unique}); do
gawk -i inplace -v name=${genome} 'BEGIN{FS=OFS="\t"} {gsub(">", ">"name";")} 1' \
${genome_dir}/${genome}.fa;
done;
done

for genome in $(cat ${genome_list}); do
for gene in $(cat ${unique}); do
gawk -i inplace -v name=${genome} 'BEGIN{FS=OFS="\t"} {gsub(">", ">"name";")} 1' \
${outdir}/${genome}/${gene}_${genome}.faa;
done;
done


#7. Concat all fastas for a particular gene across genomes-------------------------------------------------------------------------------------------------------------------------
for genome in $(cat ${genome_list}); do
for gene in $(cat ${unique}); do
mkdir -p ${outdir}/${gene}
mv ${outdir}/${genome}/${gene}_${genome}.faa ${outdir}/${gene}/;
done;
done

for genome in $(cat ${genome_list}); do
for gene in $(cat ${unique}); do
mkdir -p ${outdir}/${gene}
mv ${outdir}/${genome}/${gene}_${genome}.fna ${outdir}/${gene}/;
done;
done

for gene in $(cat ${unique}); do
cat ${outdir}/${gene}/*.faa > ${outdir}/${gene}_all_genomes.faa;
done
for gene in $(cat ${unique}); do
cat ${outdir}/${gene}/*.fna > ${outdir}/${gene}_all_genomes.fna;
done

#8. Multiple sequence alignment (written as slurm array job)-----------------------------------------------------------------------------------------------------------------------
#(written as a slurm for the Nibi Compute Canada Cluster)

mkdir -p ${outdir}/msa
module load mafft-mpi

#Nucleotide alignment
for gene in $(cat ${unique}); do
mafft --genafpair \
--maxiterate 1000 \
--thread 4 \
--reorder ${outdir}/${gene}_all_genomes.fna > ${outdir}/msa/${gene}_all_genomes_alignment.fna;
done

#Amino acid alignment
for gene in $(cat ${unique}); do
mafft --genafpair \
--maxiterate 1000 \
--thread 4 \
--reorder ${outdir}/${gene}_all_genomes.faa > ${outdir}/msa/${gene}_all_genomes_alignment.faa;
done

#9. Filter out columns with >80% gaps using TrimAl---------------------------------------------------------------------------------------------------------------------------------
#(written as a slurm array job for the Nibi Compute Canada Cluster)
module load trimal

#Filter nucleotide alignments
for gene in $(cat ${unique}); do
trimal -in ${outdir}/msa/${gene}_all_genomes_alignment.fna -out ${outdir}/msa/${gene}_all_genomes_alignment_filt.fna \
-gt 0.80 \
-sident \
-soverlap \
-htmlout ${outdir}/msa/${gene}_all_genomes_alignment_filt.html;
done

#Filter amino acid alignments
for gene in $(cat ${unique}); do
trimal -in ${outdir}/msa/${gene}_all_genomes_alignment.faa -out ${outdir}/msa/${gene}_all_genomes_alignment_filt.faa \
-gt 0.80 \
-sident \
-soverlap \
-htmlout ${outdir}/msa/${gene}_all_genomes_alignment_filt.html;
done


#10. Get counts for each gene in each genome---------------------------------------------------------------------------------------------------------------------------------------
	##We are going to use the the hmm_w_gff files we created at the beginning to get counts for each unique gene.
	##You just need to define the column that describes the gene annotation you would like to count by.

##Change NF to the column for which you want to organize/count by
mkdir -p ${outdir}/gene_counts
mkdir -p ${outdir}/gene_counts/separated_by_genome

for genome in $(cat ${genome_list}); do
awk 'BEGIN{FS=OFS="\t"} NR==1{print "gene_hit\tcount_in_genome\tgenome"; next} NR>1 {count[$(NF)]++} END {for (i in count) print i, count[i]}' ${outdir}/${genome}_hmm_w_gff.tsv > ${outdir}/gene_counts/separated_by_genome/${genome}_gene_counts.tsv

gawk -i inplace -v genome=${genome} -FS'\t' -v OFS='\t' 'NR==1{$(NF+1)="genome"; print; next} {$(NF+1)=genome; print}' \
${outdir}/gene_counts/separated_by_genome/${genome}_gene_counts.tsv;
done

##Create a summary count file (1 file per gene, but counts from genomes are combined)
	##The unique variable should reflect the exact gene annotations used to count by.
	##In my case this was a combination of the gene name and product description I made earlier on.

for genome in $(cat ${genome_list}); do
for gene in $(cat ${unique}); do
awk -F'\t' -v OFS='\t' \
    -v col=1 \
    -v val=${gene} \
    '$col==val { print }' \
${outdir}/gene_counts/separated_by_genome/${genome}_gene_counts.tsv \
>> ${outdir}/gene_counts/gene_counts_all_genomes.tsv;
done;
done

#Filtered alignments and html files are used for downstream analyses (e.g., tree construction)
