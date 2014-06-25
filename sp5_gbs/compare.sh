#!/bin/bash

#  
#  Copyright 2014 INRA-CIRAD
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or 
#  write to the Free Software Foundation, Inc., 
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

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
