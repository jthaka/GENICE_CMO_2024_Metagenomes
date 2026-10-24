#Goal: Plot the distribution of percent identity scores for hydrocarbon marker genes (R script)
#between the metagenome assembled genomes (MAGs) in our study and publicly available genomes
#belonging to the bacterial genus Cycloclasticus. 
#Required inputs: multiple sequence alignment for gene(s) of interest.

#Load required packages
if (!require('bio3d')) install.packages("bio3d"); library('bio3d')
if (!require('tidyverse')) install.packages("tidyverse"); library('tidyverse')

alignment="path/to/alignment/file.fa"

#Function for making percent identity matrix from alignment file
alignment_dist_matrix= function(fasta){
mat=bio3d::seqidentity(bio3d::read.fasta(print(fasta)))
mat[lower.tri(mat, diag = TRUE)] <- NA
mat <- mat[-nrow(mat), -1]
mat %>%
as.data.frame %>%
rownames_to_column("seq") %>%
gather(seq2, pid, - 1) %>%
mutate(seq = fct_inorder(seq),
seq2 = fct_inorder(seq2),
pid2 = ifelse(!is.na(pid), paste0(pid * 100, "%"), NA_character_)) %>%
ggplot()+
geom_tile(aes(x = seq, y = seq2, fill = pid)) +
scale_fill_gradient(low = "white",
high = "13ae4b",
na.value = "white",
breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
limits = c(0, 1)) +
geom_text(aes(x = seq, y = seq2, label = pid2)) +
theme_classic() +
scale_x_discrete(position = "top") +
theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
xlab("") +
ylab("") +
labs(title =paste("% identity of fasta between internal and external Cycloclasticus MAGs/Genomes for",substitute(fasta)))+
coord_equal()
}

#Function for making frequency plot of percent identity values from alignment file
histogram_of_alignment= function(fasta){
mat=bio3d::seqidentity(bio3d::read.fasta(print(fasta)))
mat[lower.tri(mat, diag = TRUE)] <- NA
mat <- mat[-nrow(mat), -1]
nseq=nrow(mat)
freq=table(mat)
barplot(freq, xlab="Identity between two sequences", ylab="Frequency of identity value", main=paste("Distribution of percent ID values between internal and external cycloclasticus genomes for",substitute(fasta),"(#sequences=",nseq,")"))
        }

#Make the plots and export
pdf(file = "path/to/output/for/matrix_file.pdf",
 width = 50,
height = 50)
alignment_dist_matrix(AlkB)
dev.off()

pdf(file = "path/to/output/for/histogram_file.pdf",
 width = 15,
height = 10)
histogram_of_alignment(AlkB)
dev.off()

