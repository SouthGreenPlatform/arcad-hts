#!/bin/bash

#  
#  Copyright 2014 INRA-CIRAD
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
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


for i in `ls *_reverse.fastq`;do

ref=${i%_reverse.fastq}
forward="/home/chouiki/work/Olivier/compare/"$ref"_reverse.fastq"
reverse="/home/chouiki/work/Olivier/compare/"$ref"_forward.fastq"
pairedfasta=`echo $i | cut -d_ -f4`
paired=`echo $pairedfasta | cut -d. -f1`
individu=`echo $ref | cut -d- -f2`
index=`echo $i | cut -d_ -f2`
id=$individu
tab='\t'

echo -e $forward$tab$id$tab"paired"$tab$reverse >> ../mapping/conf.txt;
done;

for i in `ls *_single.fastq`;do

ref=${i%_single.fastq}
single="/home/chouiki/work/Olivier/compare/"$ref"_single.fastq"
pairedfasta=`echo $i | cut -d_ -f4`
paired=`echo $pairedfasta | cut -d. -f1`
individu=`echo $ref | cut -d- -f2`
index=`echo $i | cut -d_ -f2`
id=$individu
tab='\t'
echo -e $single$tab$id$tab"single" >> ../mapping/conf.txt;

done;

echo "DONE.";

