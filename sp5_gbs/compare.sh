#!/bin/bash


#mkdir compare 
#CB4_c_TAGCTT_L001_R1.fastq.filtered
#OS_S1_GTTTCG_L001_-OS_015_1.fastq.filtered
#S1-Vd002_1.fastq.filtered
#S1-VC081_1.fastq.filtered


for i in `ls *1.fastq.filtered`;do
ref=${i%1.fastq.filtered}

inpath2=$ref"2.fastq.filtered"
forward="../compare/"$ref"forward.fastq"
reverse="../compare/"$ref"reverse.fastq"
single="../compare/"$ref"single.fastq"

qsub -q bioinfo.q -b yes -V -cwd -N "compare" perl /NAS/arcad_data/Softs/trunk/3_compare_fastq_paired_v5.pl -f $i -r $inpath2 -of $forward -or $reverse -os $single

done;

echo "DONE.";
