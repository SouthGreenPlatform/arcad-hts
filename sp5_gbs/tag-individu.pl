#!/usr/bin/perl -w


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

use strict;


my $i = 0;
my $j = 1;

my $output1 = "/home/mroques/DATA/Miseq1/data_dezip";
open (OUT1, $output1) or die "erreur ouverture fichier $output1";
my $output2 ="/home/mroques/DATA/Miseq1/data_dezip/S$j";
open (OUT2, ">$output2") or die "erreur ouverture fichier $output2";

while (my $lignes = <OUT1>){
		chomp $lignes;
		my @table = split (/\t/, $lignes);
		my $individu = $table[3];
		my $tag = $table[5];
		if ($i <= 7){
			print OUT2 "$tag\t$individu\n";
			$i++;
		}

		else{
			print OUT2 "*\trobus\n";
			$j++;
			my $output2 ="/home/mroques/DATA/Miseq1/data_dezip/S$j";
			open (OUT2, ">$output2") or die "erreur ouverture fichier $output2";
			print OUT2 "$tag\t$individu\n";
			$i = 1;
		}
}
	print OUT2 "*\trobus\n";

close OUT1;
close OUT2;

