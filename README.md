## Method Development: Using short read data to call TEs as present/absent

A number of different parameters were tested to see what affects TE calling accuracy. 

B73, Mo17, PH207 and W22 short reads were aligned to the B73 reference genome. Homologous TE presence/absence calls in B73, Mo17, PH207 and W22 obtained from Sarah Anderson. Dataset is described in this paper: https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.14489

The effects of (1) TE fragment length (e.g. 1bp or 3bp around the start and end of the TE), (2) the genome wide mean coverage of the sample bam file, (3) coverage before and after the start and end of a TE versus coverage just after the start and just before the end of a TE, and (4) the minimum coverage over a TE fragment neeed to call a TE as present or absent on TE calling accuracy were tested. 

bedtools multicov (v2.29.2) was used to get coverage over TE fragments. 

Here are example bedtools commands. 

```
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.20.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_1bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.20_start.in.multicov_1bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.32.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_1bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.32_start.in.multicov_1bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.45.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_1bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.45_start.in.multicov_1bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.20.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_3bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.20_start.in.multicov_3bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.32.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_3bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.32_start.in.multicov_3bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.45.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_3bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.45_start.in.multicov_3bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.20.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_5bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.20_start.in.multicov_5bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.32.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_5bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.32_start.in.multicov_5bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.45.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_5bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.45_start.in.multicov_5bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.20.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_10bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.20_start.in.multicov_10bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.32.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_10bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.32_start.in.multicov_10bp.txt
bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.45.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_10bp.gff3 > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/Mo17_refB73v4_subsample0.45_start.in.multicov_10bp.txt
```

Seperate bedtools multicov files were created for each coverage at each end of the TE (just before and after the start and end of the TE). That is, for every unique genotype, subsample and TE fragment length combination 4 bedtools files were created. 

Those files were combined and then input into R to get true and false positive/negative TE calling rates for different coverage thresholds. 
