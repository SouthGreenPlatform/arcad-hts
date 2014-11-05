#!/bin/bash

#cmp_or_quit expected returned
function cmp_or_quit {
    cmp $1 $2
    if [[ "$?" == 0 ]]; then
        rm $2
    else 
        exit
    fi
}  

echo "test paired:"
python2.6 ../demultadapt.py -l 1.0 -f indi_A_1.fastq -F indi_A_2.fastq -p returned -v adaptateur.txt
echo "check:"
cmp_or_quit expected-indiv_1_1.fastq returned-indiv_1_1.fastq 
cmp_or_quit expected-indiv_1_2.fastq returned-indiv_1_2.fastq
cmp_or_quit expected-indiv_2_1.fastq returned-indiv_2_1.fastq
cmp_or_quit expected-indiv_2_2.fastq returned-indiv_2_2.fastq 
cmp_or_quit expected-rebut_1.fastq   returned-rebut_1.fastq
cmp_or_quit expected-rebut_2.fastq   returned-rebut_2.fastq 

echo "TEST PASS"


echo "test single"
python2.6 ../demultadapt.py -l 1.0 -f indi_A_1.fastq -p returned -v adaptateur.txt
echo "check:"
cmp_or_quit expected-indiv_1.fastq returned-indiv_1.fastq 
cmp_or_quit expected-indiv_2.fastq returned-indiv_2.fastq
cmp_or_quit expected-rebut.fastq   returned-rebut.fastq

echo "TEST PASS"
