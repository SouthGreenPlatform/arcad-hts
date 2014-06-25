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

