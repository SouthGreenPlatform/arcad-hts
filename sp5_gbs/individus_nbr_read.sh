#!/bin/bash
# premier

##################################################################################################################
#Ce script parse les répértoires générer par FastqC et donne à la fin le nombre de reads existant dans fichiers
##################################################################################################################


#RPI2_S1_L001_R1_001_fastqc

for i in *1_fastqc;do 
cd $i

awk '/Total Sequences/ { printf("%s ", $3); next }' fastqc_data.txt  >> ./readsParIndividus.txt
awk '/Filename/ { print $2 }' fastqc_data.txt >> ./readsParIndividus.txt

cd ..
done;

echo "DONE.";

