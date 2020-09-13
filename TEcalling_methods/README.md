# TE_pop_gen
Scripts used to call TE insertions and analyze TE insertion patterns in the WiDiv panel. 

TEidentification folder has scripts used to generate commands to run bedtools multicov in Mo17, PH207 and W22 reads aligned to B73 over gff files that have the coordinates of different bp lengths around the start and end of every TE in the B73 genome. 

TE annotation gff3 files obtained here: https://mcstitzer.github.io/maize_TEs/
The full length gff3 files were used. 

bedtools/2.29.2 was used to obtain mean coverage

samtools/1.9 was used to subsample Mo17 and PH207 bams to lower coverage levels 

R/3.6.0 was used to call TE presence and absence

R/3.6.2 used to analyze TE presence and absence calling rates in reference genome datasets, run models on the final TE calling dataset and to make figures.

TE presence/absence calls in B73, Mo17, PH207 and W22 obtained from Sarah Anderson. Dataset is described in this paper: https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.14489

Gene expression and TE expression data also obtained from Sarah Anderson. Dataset is described in this paper: https://www.g3journal.org/content/9/11/3673.long

TE_analysis folder has scripts used to analyze TE precense/absence calls. 
