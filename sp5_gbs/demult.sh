#!/bin/bash

#RPi2_S1_L001_R1_001.fastq



for i in `ls *R1_001.fastq`;do
ref=${i%R1_001.fastq}
numIndex=`echo $i | cut -d_ -f2`
forward="/home/mroques/DATA/Miseq5/data_dezip/"$ref"R1_001.fastq"
reverse="/home/mroques/DATA/Miseq5/data_dezip/"$ref"R2_001.fastq"
#prefix=`echo $i | cut -d_ -f1-2`
fileIndex="/home/mroques/DATA/Miseq5/data_dezip/"$numIndex


qsub -q bioinfo.q -b yes -V -cwd -N "demult"$numIndex python26 /home/sarah1/Scripts/Vincent/demultadapt/demultadapt.py -f $forward -F $reverse -p $numIndex $fileIndex

done;

echo "DONE.";

