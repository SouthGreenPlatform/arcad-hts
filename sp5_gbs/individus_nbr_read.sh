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

