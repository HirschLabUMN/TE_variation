################## B73 section 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qiuxx221@umn.edu
#SBATCH -o parse_b73ref.out
#SBATCH -e parse_b73ref.err


module load bedtools 
module load samtools 

cd /scratch.global/qiuxx221/widi_bam_self/B73

for i in *.bam; do
  bedtools multicov -bams $i -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_10bp.gff3  | cut -f 9,10 | sed 's|=|\t|g;s|;|\t|g' | cut -f 2,5 | sed 's|$|\tstart|'>  $(basename "$i")_start_10bp.txt
  bedtools multicov -bams $i -bed /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_end.in_10bp.gff3  | cut -f 9,10 | sed 's|=|\t|g;s|;|\t|g' | cut -f 2,5 | sed 's|$|\tend|'>  $(basename "$i")_end_10bp.txt
  cat $(basename "$i")_start_10bp.txt $(basename "$i")_end_10bp.txt > $(basename "$i")_combined_10bp.txt
done

for j in *_combined_10bp.txt; do
    id=$(echo "$j" | cut -d'_' -f4 | cut -d '.' -f 1)
    seq_genome=$(echo "$j" |cut -d '_' -f2)
    ref_genome=$(echo "$j" |cut -d '_' -f3)
    sed "s/$/\t$id/;s/$/\t$seq_genome/;s/$/\t$ref_genome/" $j | sed "1ite_name\t10bp_coverage\tside\tgenome_cov\tseq_genome\tref_genome" > $(basename "$j")_fmt_col_name.txt
done


################# Mo17 section 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qiuxx221@umn.edu
#SBATCH -o parse_Mo17ref.out
#SBATCH -e parse_Mo17ref.err

module load bedtools 
module load samtools 
cd /scratch.global/qiuxx221/widi_bam_self/Mo17 

for i in *.bam; do
  bedtools multicov -bams $i -bed /home/hirschc1/oconnorc/sarah_te_files/Mo17.structuralTEv2.2018.12.28.filteredTE_start.in_10bp.gff3  | cut -f 9,10 | sed 's|=|\t|g;s|;|\t|g' | cut -f 2,5 | sed 's|$|\tstart|'>  $(basename "$i")_start_10bp.txt
  bedtools multicov -bams $i -bed /home/hirschc1/oconnorc/sarah_te_files/Mo17.structuralTEv2.2018.12.28.filteredTE_end.in_10bp.gff3  | cut -f 9,10 | sed 's|=|\t|g;s|;|\t|g' | cut -f 2,5 | sed 's|$|\tend|'>  $(basename "$i")_end_10bp.txt
  cat $(basename "$i")_start_10bp.txt $(basename "$i")_end_10bp.txt > $(basename "$i")_combined_10bp.txt
done

for j in *_combined_10bp.txt; do
    id=$(echo "$j" | cut -d'_' -f4 | cut -d '.' -f 1)
    seq_genome=$(echo "$j" |cut -d '_' -f2)
    ref_genome=$(echo "$j" |cut -d '_' -f3)
    sed "s/$/\t$id/;s/$/\t$seq_genome/;s/$/\t$ref_genome/" $j | sed "1ite_name\t10bp_coverage\tside\tgenome_cov\tseq_genome\tref_genome" > $(basename "$j")_fmt_col_name.txt
done





################# PH207 section 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qiuxx221@umn.edu
#SBATCH -o parse_PH207ref.out
#SBATCH -e parse_PH207ref.err

module load bedtools 
module load samtools 
cd /scratch.global/qiuxx221/widi_bam_self/PH207

for i in *.bam; do
  bedtools multicov -bams $i -bed /home/hirschc1/oconnorc/sarah_te_files/PH207.structuralTEv2.2018.10.30.filteredTE_start.in_10bp.gff3  | cut -f 9,10 | sed 's|=|\t|g;s|;|\t|g' | cut -f 2,5 | sed 's|$|\tstart|'>  $(basename "$i")_start_10bp.txt
  bedtools multicov -bams $i -bed /home/hirschc1/oconnorc/sarah_te_files/PH207.structuralTEv2.2018.10.30.filteredTE_end.in_10bp.gff3  | cut -f 9,10 | sed 's|=|\t|g;s|;|\t|g' | cut -f 2,5 | sed 's|$|\tend|'>  $(basename "$i")_end_10bp.txt
  cat $(basename "$i")_start_10bp.txt $(basename "$i")_end_10bp.txt > $(basename "$i")_combined_10bp.txt
done

for j in *_combined_10bp.txt; do
    id=$(echo "$j" | cut -d'_' -f4 | cut -d '.' -f 1)
    seq_genome=$(echo "$j" |cut -d '_' -f2)
    ref_genome=$(echo "$j" |cut -d '_' -f3)
    sed "s/$/\t$id/;s/$/\t$seq_genome/;s/$/\t$ref_genome/" $j | sed "1ite_name\t10bp_coverage\tside\tgenome_cov\tseq_genome\tref_genome" > $(basename "$j")_fmt_col_name.txt
done



################# W22 section 

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qiuxx221@umn.edu
#SBATCH -o parse_W22ref.out
#SBATCH -e parse_W22ref.err

module load bedtools 
module load samtools 
cd /scratch.global/qiuxx221/widi_bam_self/W22


for i in *.bam; do
  bedtools multicov -bams $i -bed /home/hirschc1/oconnorc/sarah_te_files/W22.structuralTEv2.10.08.2018.filteredTE_start.in_10bp.gff3  | cut -f 9,10 | sed 's|=|\t|g;s|;|\t|g' | cut -f 2,5 | sed 's|$|\tstart|'>  $(basename "$i")_start_10bp.txt
  bedtools multicov -bams $i -bed /home/hirschc1/oconnorc/sarah_te_files/W22.structuralTEv2.10.08.2018.filteredTE_end.in_10bp.gff3  | cut -f 9,10 | sed 's|=|\t|g;s|;|\t|g' | cut -f 2,5 | sed 's|$|\tend|'>  $(basename "$i")_end_10bp.txt
  cat $(basename "$i")_start_10bp.txt $(basename "$i")_end_10bp.txt > $(basename "$i")_combined_10bp.txt
done

for j in *_combined_10bp.txt; do
    id=$(echo "$j" | cut -d'_' -f4 | cut -d '.' -f 1)
    seq_genome=$(echo "$j" |cut -d '_' -f2)
    ref_genome=$(echo "$j" |cut -d '_' -f3)
    sed "s/$/\t$id/;s/$/\t$seq_genome/;s/$/\t$ref_genome/" $j | sed "1ite_name\t10bp_coverage\tside\tgenome_cov\tseq_genome\tref_genome" > $(basename "$j")_fmt_col_name.txt
done


mv *combined_10bp.txt_fmt_col_name.txt Mo17_ref_fmt_10_bp/
rename combined_10bp.txt_fmt_col_name.txt fmt_10_bp.txt *txt

mkdir PH207_ref_fmt_10_bp
mv *combined_10bp.txt_fmt_col_name.txt PH207_ref_fmt_10_bp
cd PH207_ref_fmt_10_bp
rename combined_10bp.txt_fmt_col_name.txt fmt_10_bp.txt *txt
cd ..
cp -r PH207_ref_fmt_10_bp ~/TE_variation



mkdir W22_ref_fmt_10_bp
mv *combined_10bp.txt_fmt_col_name.txt W22_ref_fmt_10_bp
cd W22_ref_fmt_10_bp
rename combined_10bp.txt_fmt_col_name.txt fmt_10_bp.txt *txt
cd ..
cp -r W22_ref_fmt_10_bp ~/TE_variation




