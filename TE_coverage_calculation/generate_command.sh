#Script to generate commands to run bedtools multicov on reference genome test samples 

ref_names="B73
Mo17
PH207
W22"                                                                                                                                                                                                                                                             

cov_len="1
3
5
10"

#generate bedtools multicov command to get coverage of TE start and end fragments of different length
for i in $ref_names
do
  for num in $cov_len
    do
      echo "bedtools multicov -bams /panfs/roc/scratch/oconnorc/reference_bams/${i}_B73v4_map.q20.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.out_${num}bp.gff3 -s > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_start.out.multicov_${num}bp.txt" >> $1
      echo "bedtools multicov -bams /panfs/roc/scratch/oconnorc/reference_bams/${i}_B73v4_map.q20.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_${num}bp.gff3 -s > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_start.in.multicov_${num}bp.txt" >> $1
      echo "bedtools multicov -bams /panfs/roc/scratch/oconnorc/reference_bams/${i}_B73v4_map.q20.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_end.out_${num}bp.gff3 -s > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_end.out.multicov_${num}bp.txt" >> $1
      echo "bedtools multicov -bams /panfs/roc/scratch/oconnorc/reference_bams/${i}_B73v4_map.q20.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_end.in_${num}bp.gff3 -s > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_end.in.multicov_${num}bp.txt" >> $1
    done
done

subsample_ref="B73
Mo17"
subsample_lvl="0.20
0.32
0.45
0.66
0.83"

#generate bedtools multicov command to get coverage of TE start and end fragments of different length for Mo17 and PH207 bam files subsampled to lower coverage levels
for i in $subsample_ref
  do
    for num in $cov_len
      do
        for lvl in $subsample_lvl
          do
            echo "bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled/${i}_B73v4_map.q20_subsample${lvl}.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.out_${num}bp.gff3 -s > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_start.out.multicov_${num}bp_subsample${lvl}.txt" >> $1
            echo "bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled/${i}_B73v4_map.q20_subsample${lvl}.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_${num}bp.gff3 -s > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_start.in.multicov_${num}bp_subsample${lvl}.txt" >> $1
            echo "bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled/${i}_B73v4_map.q20_subsample${lvl}.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_end.out_${num}bp.gff3 -s > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_end.out.multicov_${num}bp_subsample${lvl}.txt" >> $1
            echo "bedtools multicov -bams /home/hirschc1/oconnorc/TE_project/subsampled/${i}_B73v4_map.q20_subsample${lvl}.bam -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_end.in_${num}bp.gff3 -s > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_end.in.multicov_${num}bp_subsample${lvl}.txt" >> $1
          done
      done
  done
#What are sub-sampled files? Mo17 and PH207 both have genome wide mean coverage of ~30x but W22 has a genome wide mean coverage of 6x and most samples have genome wide mean coverage between 10x - 20x.
#I wanted to look at TE calling rate at different coverage levels, so I used samtools view to create lower coverage versions of the same bam files
#command used to sub-sample: samtools view -s 0.66 /panfs/roc/scratch/oconnorc/reference_bams/Mo17_B73v4_map.q20.bam -O bam -o /home/hirschc1/oconnorc/TE_project/subsampled_bams/Mo17_B73v4_map.q20_subsample0.66.bam

#example of a command used to join the 4 TE sides together: 
#echo \"TE_name    TE_status       start.out       start.in        end.in  end.out\" > /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_TEcov_${num}bp_${i}status.txt" >> $1
#join -1 1 -2 1 <(cut -f 9,10 /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_start.out.multicov_${num}bp.txt | sort -k 1,1) <(cut -f 9,10 /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_start.in.multicov_${num}bp.txt | sort -k 1,1) | tr \" \" \"\t\" | join -1 1 -2 1 - <(cut -f 9,10 /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_end.in.multicov_${num}bp.txt | sort -k 1,1) | tr \" \" \"\t\" | join -1 1 -2 1 - <(cut -f 9,10 /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_end.out.multicov_${num}bp.txt | sort -k 1,1) | tr \" \" \"\t\" | tr \"=\" \"\t\" | tr \";\" \"\t\" | cut -f 2,5-8 | join -1 1 -2 1 <(sort -k 1,1 /home/hirschc1/oconnorc/TE_project/B73_TEs_${i}status.txt) <(sort -k 1,1 -) | tr \" \" \"\t\" >> /home/hirschc1/oconnorc/TE_project/bedtools_multicov_test/${i}_refB73v4_TEcov_${num}bp_${i}status.txt

#end of commands 
