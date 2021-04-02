# This study used random forest machine learning models and 509 whole genome resequencing data to predict the presence/absence of a known TE in a diverse maize panel. 

TEidentification folder has scripts used to generate commands to run bedtools multicov in Mo17, PH207 and W22 reads aligned to B73 over gff files that have the coordinates of different bp lengths around the start and end of every TE in the B73 genome. 

TE annotation gff3 files obtained here: https://mcstitzer.github.io/maize_TEs/ The full length gff3 files were used. 

bedtools/2.29.2 was used to obtain mean coverage
sambamba 0.8.0 was used to downsample bams to lower coverage levels (15x, and 30x) for machine learning traning model. If coverage was at or below either of these thresholds no down sampling was performed. 

TE presence/absence calls in B73, Mo17, PH207 and W22 obtained from Sarah Anderson. Dataset is described in this paper: https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.14489

The training set consisted of counts from resequencing data from three of the four genomes aligned to all four reference genome assemblies and the test set consisted of the counts from resequencing data from the fourth genome mapped to all for reference genome assemblies.The mean coverage over the 10 bp at the start of the TE, mean coverage over the 10 bp at the end of the TE, and the TE order were used as predictors to train the model. Separate models were trained for realized coverage of 30x and realized coverage of 15x. 

For widiv panel, if the realized coverage for a sample was >=25x depth the 30x model was used and if the realized coverage for a sample was <25x depth the 15x model was used. If the probability of present from the model was >=0.7 a TE was classified as present in the sample. If the probability of present from the model was <=0.3 a TE was classified as absent in the sample. All other TEs were classified as ambiguous. 

For any TE where the resequencing reads mapped to its cognate genome (e.g. B73 reads mapped to the B73 genome assembly) did not result in a present classification the TE was considered recalcitrant to accurate calls for short read data and was removed from downstream analysis. Across all samples mapped to a reference genome assembly if there was greater than 25% ambiguous data for a TE the TE was removed from downstream.

A final non-redundant TE matrix was merged using homology TE matrix from the (Anderson et al., 2019)

Fhe final TE varation matrix for the Widiv panel includes TE IDs, TE population frequency, TE Family size, LTR similarity, and nested status. 

R/3.6.3 used to analyze TE presence and absence calling and downstream figure visulization. 




