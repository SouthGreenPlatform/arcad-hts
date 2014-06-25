#!/bin/env perl

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
=pod

=head1 NAME

vcfTOmarker.pl - transforming 01 genotype coding (vcf file) into nucleotide genotype coding (csv file) and assemble phased SNPs

=head1 DESCRIPTION

Transforming 01 genotype coding into nucleotide genotype coding
 
and assemble phased SNPs

If a genotypic class is present only once, transform it into missing data

=head1 SYNOPSIS / USAGE

vcfTOmarker.pl	[--help] [--man] [--version]
[-i VCF-file] \
[-o output-file] \
[-L fragment size] \

=cut

use strict;
use Pod::Usage;
use Getopt::Long;
use Carp qw (cluck confess croak);


=head1 OPTIONS

=over 8

=item B<--help>
Show the brief help information.

=item B<--man>
Read the manual, with examples.

=item B<--version>
Show the version number and exit.

=head2 required parameters

=over 5

=item B<[-i]> ([input_file]): 

VCF file


=item B<[-o]> ([output_file]): 

output file

=item B<[-L]> ([size]): 

fragment size

=cut

sub help
{ 
	my $msg = shift;
	
	print "-----------------------------------------------------------------------------\n"
	."-- ERROR MESSAGE --\n"."$msg\n".
	"-----------------------------------------------------------------------------\n"
	if( defined  $msg );
	pod2usage(0);
}

help "Not enough arguments" unless (@ARGV);

my ($man, $help, $version, $input, $output, $taille);
$man = 0;
$help = 0;
$version = 0;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"input|i=s"   => \$input,
	"output|o=s"  => \$output,
	"taille|L=s"  => \$taille	
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
if ($version) {print "1.0.0\n"; exit;}

my $output_TMP = "./TMP.csv";

open (F, $input) or die "cannot open the file $input\n";
open(OUT, ">$output_TMP") or die"cannot open the file $output_TMP\n";
open(OUT1, ">$output") or die"cannot open the file $output\n";

my (@tab, %count, @tab_base_alt, @distf, %CLASSES, %HOMO, %HETERO , %hash, %hash2, @chr, @al1, @pos, @AL1, @AL2, $dist, $GN, $phas, $nbr_mm_class, $chr, $pos, $base_ref, $base_alt, $allele1, $allele2, $alt_base1, $alt_base2, $alt_base3, $b1, $b2, $nophas, $NN, $countHOMO, $countHETERO, %nbr_HETERO, %nbr_HOMO, $classe, $nbr_classes, $total_indiv, $nbr_hetero, $nbr_homo);

$nophas = 0;
$phas = 0;
$NN =0;
$nbr_classes =0;
$nbr_mm_class =0;
$dist = 0;
%count = ( );
%hash = ();
%hash2 = ();
my $j = 0;
#Entete
print OUT "Chromosome	POS	";
while(my $line =<F>){
	chomp $line;
	if($line =~ m/^#/){
		if ($line =~ m/^#CHROM/) {
			@tab = split (/\t/, $line);
			#Suite de l'entete
			for(my $i = 9;$i<@tab;$i++){
				print OUT "$tab[$i]\t";
			}
			print OUT "\n";
		}
	}
	else{
		@tab = split (/\t/, $line);
		push (@chr, $tab[0]);
		push (@pos, $tab[1]);
		$base_ref = $tab[3];
		$base_alt = $tab[4];
		@tab_base_alt = split (',', $base_alt);
		$alt_base1 = $tab_base_alt[0];
		$alt_base2 = $tab_base_alt[1];
		$alt_base3 = $tab_base_alt[2];
		# garder que les SNPs bialléliques
		if(@tab_base_alt == 1){
			#meme chromosome
			if($chr[$j] eq $chr[$j-1]){
				$dist = $pos[$j] - $pos[$j-1];
				#print TEST "$pos[$j]**$pos[$j-1]==$dist\n";
				for(my $i = 9;$i<@tab;$i++){
				#Vérifier si tous les individus sont phasés pour un SNP donné
					if (($tab[$i] =~ m/^(\d)\/(\d)*/)){
						
					}
					elsif (($tab[$i] =~ m/^(\d)\|(\d)*/) or ($tab[$i] =~ m/^(\.)\/(\.)/)){
						$phas++;
					}
				}
				# Si le SNP est phasé le mettre dans un hash
				if($phas == ("$#tab" - 8)){
					#print OUT "$pos[$j]\n";
					push(@distf, $pos[$j]);
					my $t = $distf[-1]-$distf[0];
					#print DIST "$t\n";
					if ($t < $taille){
						for(my $i = 9;$i<@tab;$i++){
							if ($tab[$i] =~ m/^(\d)\|(\d)*/){
								my $b1;
								my $b2;
								if($1 == 0){
									$b1 = $base_ref;
								}
								elsif($1 == 1){
									$b1 = $alt_base1;
								}
								if($2 == 0){
									$b2 = $base_ref;
								}
								elsif($2 == 1){
									$b2 = $alt_base1;
								}
								# verifier si la distance du fragment est inférieiur à la taille donnée en parametre
								push(@{$hash{$i}},$b1);
								push(@{$hash2{$i}},$b2);
							}
							elsif ($tab[$i] =~ m/^(\.)\/(\.)/){
								$b1 = "N";
								$b2 = "N";
								push(@{$hash{$i}},$b1);
								push(@{$hash2{$i}},$b2);
							}
						}
					}
					#Un autre fragment 
					else{
						if(%hash eq 0){
								
						}
						else{
							print OUT "$chr[$j-1]	$pos[$j-1]\t";
						}
						foreach my $i(sort {$a<=>$b} keys (%hash)) {
							my $marker = "@{$hash{$i}}|@{$hash2{$i}}";
							$marker =~ s/ //gm;
							if($marker =~ m/N/){
								$marker = "N|N";
							}
							print OUT "$marker\t";
						}
						print OUT "\n";
						%hash = ();
						%hash2 = ();
						@distf =  ();
						push(@distf, $pos[$j]);
						for(my $i = 9;$i<@tab;$i++){
							if (($tab[$i] =~ m/^(\d)\|(\d)*/) or ($tab[$i] =~ m/^(\d)\/(\d)*/)){
								my $b1;
								my $b2;
								if($1 == 0){
									$b1 = $base_ref;
								}
								elsif($1 == 1){
									$b1 = $alt_base1;
								}
								if($2 == 0){
									$b2 = $base_ref;
								}
								elsif($2 == 1){
									$b2 = $alt_base1;
								}
								#Conserver la ligne
								push(@{$hash{$i}},$b1);
								push(@{$hash2{$i}},$b2);
							}
							elsif ($tab[$i] =~ m/^(\.)\/(\.)/){
								$b1 = "N";
								$b2 = "N";
								push(@{$hash{$i}},$b1);
								push(@{$hash2{$i}},$b2);
							}
						}
					}
				}
				# Pas phasé
				else{
					if(%hash eq 0){
							
					}
					else{
						print OUT "$chr[$j-1]	$pos[$j-1]\t";
					}
					foreach my $i(sort {$a<=>$b} keys (%hash)) {
						my $marker = "@{$hash{$i}}|@{$hash2{$i}}";
						$marker =~ s/ //gm;
						if($marker =~ m/N/){
							$marker = "N|N";
						}
						print OUT "$marker\t";
					}
					print OUT "\n";
					%hash = ();
					%hash2 = ();
					@distf =  ();
					push(@distf, $pos[$j]);
					for(my $i = 9;$i<@tab;$i++){
						if (($tab[$i] =~ m/^(\d)\|(\d)*/) or ($tab[$i] =~ m/^(\d)\/(\d)*/)){
							my $b1;
							my $b2;
							if($1 == 0){
								$b1 = $base_ref;
							}
							elsif($1 == 1){
								$b1 = $alt_base1;
							}
							if($2 == 0){
								$b2 = $base_ref;
							}
							elsif($2 == 1){
								$b2 = $alt_base1;
							}
							#Conserver la ligne
							push(@{$hash{$i}},$b1);
							push(@{$hash2{$i}},$b2);
						}
						elsif ($tab[$i] =~ m/^(\.)\/(\.)/){
							$b1 = "N";
							$b2 = "N";
							push(@{$hash{$i}},$b1);
							push(@{$hash2{$i}},$b2);
						}
					}
				}
			}
			# Chromosome different
			else{
				#écrire ce qui a dans le hash et pushé le nouveau SNP pour verifié si le suivant est phasé ou pas
				if(%hash eq 0){
						
				}
				else{
					print OUT "$chr[$j-1]	$pos[$j-1]\t";
				}
				foreach my $i(sort {$a<=>$b} keys (%hash)) {
					my $marker = "@{$hash{$i}}|@{$hash2{$i}}";
					$marker =~ s/ //gm;
					if($marker =~ m/N/){
						$marker = "N|N";
					}
					print OUT "$marker\t";
				}
				print OUT "\n";
				%hash = ();
				%hash2 = ();
				@distf =  ();
				push(@distf, $pos[$j]);
				for(my $i = 9;$i<@tab;$i++){
					if (($tab[$i] =~ m/^(\d)\|(\d)*/) or ($tab[$i] =~ m/^(\d)\/(\d)*/)){
						my $b1;
						my $b2;
						if($1 == 0){
							$b1 = $base_ref;
						}
						elsif($1 == 1){
							$b1 = $alt_base1;
						}
						if($2 == 0){
							$b2 = $base_ref;
						}
						elsif($2 == 1){
							$b2 = $alt_base1;
						}
						#Conserver la ligne
						push(@{$hash{$i}},$b1);
						push(@{$hash2{$i}},$b2);
					}
					elsif ($tab[$i] =~ m/^(\.)\/(\.)/){
						$b1 = "N";
						$b2 = "N";
						push(@{$hash{$i}},$b1);
						push(@{$hash2{$i}},$b2);
					}
				}
			}
		}
		$j++;
		$phas = 0;
	}
}
print OUT "$chr[$j-1]	$pos[$j-1]\t";

foreach my $i(sort {$a<=>$b} keys (%hash)) {
	my $marker = "@{$hash{$i}}|@{$hash2{$i}}";
	$marker =~ s/ //gm;
	if($marker =~ m/N/){
		$marker = "N|N";
	}
	print OUT "$marker\t";
}
print OUT "\n";

close F;
close OUT;
#close TEST;


open(OUT, $output_TMP) or die"cannot open the file $output_TMP\n";

$nbr_classes = 0;
$nbr_homo = 0;
$nbr_hetero = 0;
$countHOMO = 0;
$countHETERO = 0;
$NN = 0;
my $NA = 0;
%count = ();
#my $head =<OUT>;
#print OUT1 "$head";

while(my $line =<OUT>){
	chomp $line;
	my @tab = split (/\t/, $line);
	if($line =~ m/^Chromosome/){
		@tab = split (/\t/, $line);
		for(my $i = 0;$i<@tab;$i++){
			print OUT1 "$tab[$i]\t";
		}
		print OUT1 "\n";
	}
	else{
	#calcul du nombre de classes génotypique et de l'effectif de chaque classe
	for(my $i = 4;$i<@tab;$i++){
		if($tab[$i] =~ m/^(\w*)\|(\w*)/){
			$classe = "$1|$2";
			#$classe = "$1$2";
			#print OUT1 "$classe\t";
			if($classe ne "NN"){
				if((!exists($CLASSES{"$1|$2"})) and (!exists($CLASSES{"$2|$1"}))){
					if($1 eq $2){
						$CLASSES{"$1|$2"} = 1;
					}
					else{
						$CLASSES{"$1|$2"} = 1;
						$CLASSES{"$2|$1"} = 1;
					}
					$nbr_classes++;
					if($1 eq $2){
						$HOMO{$classe} = 1;
						$nbr_homo++;
					}
					else{
						$HETERO{$classe} = 1;
						$nbr_hetero++;
					}
				}
				if((exists($CLASSES{"$1|$2"})) or (exists($CLASSES{"$2|$1"}))){
					if($1 eq $2){
						$count{"$1|$2"}++;
					}else{
						$count{"$1|$2"}++;
						$count{"$2|$1"}++;
					}
					if($1 eq $2){
						$countHOMO++;
					}
					else{
						$countHETERO++;
					}
				}
			}
		}
	}
	if ($nbr_classes > 1){
		print OUT1 "$tab[0]\t$tab[1]\t$tab[2]\t$tab[3]\t";
		for(my $i = 4;$i<@tab;$i++){
			if ($tab[$i] =~ m/^(\w*)\|(\w*)/){
				$classe = "$1|$2";
				# Si une classe génotypique se trouve une seule fois, la transformer en NN 
				if((($count{"$1|$2"} == 1) or ($count{"$2|$1"} == 1)) and ($classe ne "N|N")){
					$classe = "N|N"; 
					$NN++;
					#print "$line\n";
				}
				print OUT1 "$classe\t";
			}
		}
		print OUT1 "\n";
		if($NN >= 1){
			$NA++;
		}
		$NN = 0;
	}
	else{
		#print "$nbr_classes	$line\n";
	}
	#print OUT1 "$nbr_classes\n";
	$nbr_classes = 0;
	%CLASSES = ();
	%HOMO = ();
	%HETERO = ();
	$countHOMO = 0;
	$countHETERO = 0;
	%count = ();

	$nbr_hetero = 0;
	$nbr_homo = 0;
	}
}
#print "NN : $NA\n";

close OUT;
close OUT1;



=head1 INCOMPATIBILITIES

Fully compatble with any perl version

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 2

=item Hajar CHOUIKI (INRA), hajar.chouiki@supagro.inra.fr

=back

=head1 VERSION

1.0.0

=head1 DATE

25.06.2014

=head1 LICENSE AND COPYRIGHT

  Copyright 2014 INRA-CIRAD
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, if not, see <http://www.gnu.org/licenses/> 
  or write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston,
  MA 02110-1301, USA.

=cut
