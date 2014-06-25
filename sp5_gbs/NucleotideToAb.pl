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

NucleotideToAb.pl - transforms nucleotidic coding into abcd coding 

=head1 DESCRIPTION

Transforms nucleotidic coding into abcd coding 

Takes into account the genotyping errors and the case of null allele

Generates a file with eliminated SNPs 

=head1 SYNOPSIS / USAGE

NucleotideToAb.pl	[--help] [--man] [--version]

[-i tabulated-file] \
[-o output_file] \

example of input file :

Chromosome	POS	REF	ALT	female	male	Vd011	Vd012	Vd013	Vd015	...

chr1	192497	G	A	AG	GG	GG	GG	GG	AG	...

chr1	235688	G	A	AG	GG	GG	NN	GG	AG	...
...


=cut

use strict;
use Pod::Usage;
use Getopt::Long;
use Carp qw (cluck confess croak);


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

=pod

=head1 OPTIONS

=over 8

=item B<--help>
Show the brief help information.

=item B<--man>
Read the manual, with examples.

=item B<--version>
Show the version number and exit.

=head2 required parameters

=over 1

=item B<[-i]> ([input_file]): 
input file

=item B<[-o]> ([output_file]): 
output_file

=back

=cut

my ($man, $help, $version, $input, $output);
$man = 0;
$help = 0;
$version = 0;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"input|i=s"   => \$input,
	"output|o=s"  => \$output	
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
if ($version) {print "1.0.0\n"; exit}

if(defined($input) && ! -e $input )
{
	warn "$input not exists ";
	exit(0);
	
}

if(defined($output) && -e $output )
{
	warn "$output already exists ";
	exit(0);
	
}

open (F, $input)or die "cannot open the file $input\n";
open(OUT, ">$output") or die"cannot open the file $output\n";
my $sup = $output."_removed.txt";
open(SUP, ">$sup") or die"cannot open the file $sup\n";

my (@tab, @not, @NotOk, @clasHOMO, %HOMO_count, $contig, $pos, $male, $femelle, $maleb1, $maleb2, $femelleb1, $femelleb2, $indvb1, $indvb2 );
my (%CLASSES, %HOMO, %HETERO, $marker, $countHOMO, $countHETERO, %nbr_HETERO, %nbr_HOMO, $classe, $nbr_classes, $total_indiv, $nbr_hetero, $nbr_homo );

$total_indiv = 0;
$nbr_classes = 0;
$nbr_homo = 0;
$nbr_hetero = 0;
$countHOMO = 0;
$countHETERO = 0;

my %count = ( );
@clasHOMO = ();
%HOMO_count  = ();

print OUT "Chromosome	POS	Marker	Segregation	";

while(my $line =<F>){
	chomp $line;
	@tab = split (/\t/, $line);
	$contig = $tab[0];
	$pos = $tab[1];
	$femelle = $tab[4];
	$male = $tab[5];
	$male =~ m/(\w)(\w)/;
	$maleb1 = $1;
	$maleb2 = $2;
	$femelle =~ m/(\w)(\w)/;
	$femelleb1 = $1;
	$femelleb2 = $2;
	# Print Individuals 
	if($line =~ m/^Chromosome/){
		@tab = split (/\t/, $line);
		for(my $i = 6;$i<@tab;$i++){
			print OUT "$tab[$i]\t";
		}
		print OUT "\n";
	}
	else{
	#calculating the number of genotypic classes and the number of each class
	for(my $i = 6;$i<@tab;$i++){
		if($tab[$i] =~ m/(\w)(\w)/){
			$classe = "$1$2";
			if($classe ne "NN"){
				if(!exists($CLASSES{$classe})){
					$CLASSES{$classe} = 1;
					$nbr_classes++;
					if($1 eq $2){
						push(@clasHOMO, $classe);
						$HOMO{$classe} = 1;
						$nbr_homo++;
						$HOMO_count{$classe}++;
					}
					else{
						$HETERO{$classe} = 1;
						$nbr_hetero++;
					}
				}
				if(exists($CLASSES{$classe})){
					$count{$classe}++;
					if($1 eq $2){
						$HOMO_count{$classe}++;
						$countHOMO++;
					}
					else{
						$countHETERO++;
					}
				}
			}
		}
	}

	# retrieve the SNPs that have more than one genotypic class in the descent
	if($nbr_classes > 1){
		# case : ab x ab
		if(($femelle eq $male) and ($femelleb1 ne $femelleb2)){
			if($nbr_classes == 3){
				#print OUT "$line\n";
				print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
				for(my $i = 6;$i<@tab;$i++){
					$tab[$i] =~ m/(\w)(\w)/;
					$indvb1 = $1;
					$indvb2 = $2;
					$classe = $tab[$i];
					
					if($indvb1 ne $indvb2){
						print OUT "ab\t";
					}
					elsif($indvb1 eq $femelleb1){
						print OUT "aa\t";
					}
					elsif($indvb1 eq $femelleb2){
						print OUT "bb\t";
					}
					elsif($classe eq "NN"){
						print OUT "--\t";
					}
					else{
						print OUT "NO\t";
					}
				}
				print OUT "\n";
			}
			# we must not come across this case => Print in standard output
			elsif($nbr_classes == 2){
				print SUP "abxab2cla	$line\n";
			}
		}
		# case : aa x ab or ab x aa
		elsif($femelle ne $male){
			if ($nbr_classes == 2){
				# case : aa x ab
				# obs = real : AA x AC 
				if(($femelle ne "NN") and ($femelleb1 eq $femelleb2) and ($maleb1 ne $maleb2) and (($femelleb1 eq $maleb1) or ($femelleb1 eq $maleb2))){
					##print OUT "$line\n";
					print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
					for(my $i = 6;$i<@tab;$i++){
						$tab[$i] =~ m/(\w)(\w)/;
						$indvb1 = $1;
						$indvb2 = $2;
						$classe = $tab[$i];
						if((($classe eq $femelle) or ($indvb1 eq $indvb2)) and ($classe ne "NN")){
							print OUT "aa\t";
						}
						elsif(($classe eq $male) and ($indvb1 ne $indvb2)){
							print OUT "ab\t";
						}
						elsif($classe eq "NN"){
							print OUT "--\t";
						}
						# Just to check if there is a genotype that is wrong, it is noted "NO"
						else{
							print OUT "NO\t";
						}
					}
						print OUT "\n";
				}
				# obs : AA x CC ou (CC x AA) ,real : null allele: A0 x CC ou (null allele: CC x A0 ou err: AC x AA)) => double : 
				# We double the case when there are several real cases possible for an observed case
				# err for error
				# dm for missing data
				elsif(($femelle ne "NN") and ($male ne "NN") and ($femelleb1 eq $femelleb2) and ($maleb1 eq $maleb2)){
					#print OUT "$line\n";
					print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
					for(my $i = 6;$i<@tab;$i++){
						$tab[$i] =~ m/(\w)(\w)/;
						$indvb1 = $1;
						$indvb2 = $2;
						$classe = $tab[$i];
						if(($classe ne "NN") and ($indvb1 eq $indvb2)){
							print OUT "ab\t";
						}
						elsif($indvb1 ne $indvb2){
							print OUT "aa\t";
						}
						elsif($classe eq "NN"){
							print OUT "--\t";
						}
						else {
							print OUT "NO\t";
						}
					}
					print OUT "\n";
					##print OUT "$line\n";
					print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxaa>\t";
					for(my $i = 6;$i<@tab;$i++){
						$tab[$i] =~ m/(\w)(\w)/;
						$indvb1 = $1;
						$indvb2 = $2;
						$classe = $tab[$i];
						if(($classe ne "NN") and ($indvb1 eq $indvb2)){
							print OUT "ab\t";
						}
						elsif($indvb1 ne $indvb2){
							print OUT "aa\t";
						}
						elsif($classe eq "NN"){
							print OUT "--\t";
						}
						else {
							print OUT "NO\t";
						}
					}
					print OUT "\n";
				}
				# obs: AA x NN
				elsif(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and ($male eq "NN")){
					if($nbr_hetero == 1){
					#obs: AA x NN, real :(dm: AA x AC ou dm null allele: CC x A0 ou A0 x CC ou err+dm: AC x AA) => double
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $femelle){
								print OUT "aa\t";
							}
							elsif($indvb1 ne $indvb2){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxaa>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $femelle){
								print OUT "aa\t";
							}
							elsif($indvb1 ne $indvb2){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
					}
					# AA x NN ,real: err null allele:AC x 00
					elsif($nbr_hetero == 0){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $femelle){
								print OUT "aa\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
					}
				}
				# obs: NN x AC
				elsif(($femelle eq "NN") and ($maleb1 ne $maleb2)){
					# obs: NN x AC real:dm:  AA x AC
					if ($nbr_hetero == 1){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $male){
								print OUT "ab\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "NN") and (($indvb1 eq $maleb1) or($indvb1 eq $maleb2))){
								print OUT "aa\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
					}
					# obs: NN X AC ,real:null allele :00 x AC
					elsif($nbr_hetero == 0){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if(($classe ne "NN") and ($indvb1 eq $indvb2) and ($indvb1 eq $maleb1)){
								print OUT "ab\t";
							}
							elsif (($classe ne "NN") and ($indvb1 eq $indvb2) and ($indvb1 eq $maleb2)){ 
								print OUT "aa\t"
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
					}
					else{
						print SUP "aaxab:hetero#0/#1$line\n";
					}
				}
				# case : ab x aa ou aa x ab
				# obs = real : AC x AA 
				elsif(($male ne "NN") and ($femelle ne "NN") and ($maleb1 eq $maleb2) and ($femelleb1 ne $femelleb2) and (($maleb1 eq $femelleb1) or ($maleb1 eq $femelleb2))){
					#print OUT "$line\n";
					print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
					for(my $i = 6;$i<@tab;$i++){
						$tab[$i] =~ m/(\w)(\w)/;
						$indvb1 = $1;
						$indvb2 = $2;
						$classe = $tab[$i];
						if((($classe eq $male) or ($indvb1 eq $indvb2)) and ($classe ne "NN" )){
							print OUT "aa\t";
						}
						elsif($classe eq $femelle){
							print OUT "ab\t";
						}
						elsif($classe eq "NN"){
							print OUT "--\t";
						}
						else{
							print OUT "NO\t";
						}
					}
						print OUT "\n";
				}
				# NN x AA 
				elsif(($femelle eq "NN") and ($maleb1 eq $maleb2) and ($male ne "NN")){
					if($nbr_hetero == 1){
					# obs: NN x AA (real dm: AC x AA ou dm:null allele: CC x A0 ou A0 x CC ou err+dm: AA x AC) => double
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $male){
								print OUT "aa\t";
							}
							elsif($indvb1 ne $indvb2){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<aaxab>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $male){
								print OUT "aa\t";
							}
							elsif($indvb1 ne $indvb2){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
					}
					# NN x AA ,real: err null allele:  00 x AC
					elsif($nbr_hetero == 0){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $male){
								print OUT "aa\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne $male) and ($classe ne "NN")){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
					}
				}
				# AC x NN 
				elsif(($femelleb1 ne $femelleb2) and ($male eq "NN")){
					#obs: AC x NN real:dm: AC x AA
					if($nbr_hetero == 1){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $femelle){
								print OUT "ab\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "NN")){
								print OUT "aa\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
							print OUT "\n";
					}
					#obs: AC X NN ,real:null allele: AC x 00
					elsif($nbr_hetero == 0){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if(($classe ne "NN") and ($indvb1 eq $indvb2) and ($indvb1 eq $femelleb1)){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							elsif (($classe ne "NN") and ($indvb1 eq $indvb2) and ($indvb1 eq $femelleb2)){ 
								print OUT "aa\t"
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
					}
				}
			}
			# case : abxab or a0xab or abxa0
			elsif($nbr_classes == 3){
				if($nbr_hetero == 1){
					# case : abxab
					if($countHETERO >= 1.5 * ($countHOMO/2)){
						#obs: NN x AC real:dm:  AC x AC
						if(($femelle eq "NN") and ($maleb1 ne $maleb2)){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								if($classe eq $male){
									print OUT "ab\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne "NN") and ($indvb1 eq "$maleb1")){
									print OUT "aa\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne "NN") and ($indvb1 eq "$maleb2")){
									print OUT "bb\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "\n";
						}
						#obs: AC x NN real:dm: AC x AC
						elsif(($femelleb1 ne $femelleb2) and ($male eq "NN")){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								if($classe eq $femelle){
									print OUT "ab\t";
								}
								elsif(($indvb1 eq $indvb2) and ($indvb1 eq "$femelleb1")){
									print OUT "aa\t";
								}
								elsif(($indvb1 eq $indvb2) and ($indvb1 eq "$femelleb2")){
									print OUT "bb\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "\n";
						}
						#obs: AC x AA real: AC x AC
						elsif(($femelleb1 ne $femelleb2) and ($maleb1 eq $maleb2) and ($male ne "NN")){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								if($classe eq $femelle){
									print OUT "ab\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe eq $male) and ($classe ne "NN")){
									print OUT "aa\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $male) and ($classe ne "NN")){
									print OUT "bb\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "\n";
						}
						#obs: AA x AC real: AC x AC
						elsif(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and ($maleb1 ne $maleb2)){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								if($classe eq $male){
									print OUT "ab\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe eq $femelle) and ($classe ne "NN")){
									print OUT "aa\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
									print OUT "bb\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "\n";
						}
						else{
							print SUP "3clasParentAAxCC	$line\n";
						}
					}
					# case : a0xab or abxa0
					elsif ($countHETERO <= ($countHOMO/2)){
						# obs : AA X AC ou (AA x NN), real: null allele: A0 x AC ou (dm null allele: A0 x AC ou err+dm null allele: AC x A0) => double
						if(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and (($male eq "NN") or (($maleb1 ne $maleb2) and (($femelleb1 eq $maleb1) or ($femelleb1 eq $maleb2))))){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<a0xab>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								if($classe eq $femelle){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
									print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($femelleb1 eq $indvb1) or ($femelleb1 eq $indvb2))){
									print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
							print OUT "\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxa0>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								if($classe eq $femelle){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
									print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($femelleb1 eq $indvb1) or ($femelleb1 eq $indvb2))){
									print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "\n";
						}
						# obs : AC X AA ou (NN x AA), real: null allele: AC x A0 ou (dm null allele: AC x A0 ou err+dm null allele: A0 x AC) => double
						elsif(($maleb1 eq $maleb2) and ($male ne "NN") and (($femelle eq "NN") or (($femelleb1 ne $femelleb2) and (($maleb1 eq $femelleb1) or ($maleb1 eq $femelleb2))))){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxa0>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								if($classe eq $male){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $male) and ($classe ne "NN")){
									print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($maleb1 eq $indvb1) or ($maleb1 eq $indvb2))){
									print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
							print OUT "\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<a0xab>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								if($classe eq $male){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $male) and ($classe ne "NN")){
									print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($maleb1 eq $indvb1) or ($maleb1 eq $indvb2))){
									print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "\n";
						}
						#obs: AA x CC ou CC x AA, real: err: A0 x AC ou AC x A0 => double
						elsif(($femelleb1 eq $femelleb2) and ($maleb1 eq $maleb2) and ($male ne "NN") and ($femelle ne "NN")){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<a0xab>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								
								my $homo_1 = $clasHOMO[0];
								my $homo_2 = $clasHOMO[1];
								
								my ($homoPlus, $homoMoins);
								
								if(($indvb1 eq $indvb2) and ($classe ne "NN") and ($HOMO_count{$homo_1} > $HOMO_count{$homo_2})){
									$homoPlus = $homo_1;
								}
								else{
									$homoPlus = $homo_2;
								}
								if($classe eq $homoPlus){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne "NN") and ($classe ne $homoPlus)){
									print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($classe eq "$femelleb1$maleb1") or ($classe eq "$maleb1$femelleb1"))){
									print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
							print OUT "\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxa0>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								
								my $homo_1 = $clasHOMO[0];
								my $homo_2 = $clasHOMO[1];
								
								my ($homoPlus, $homoMoins);
								
								if(($indvb1 eq $indvb2) and ($classe ne "NN") and ($HOMO_count{$homo_1} > $HOMO_count{$homo_2})){
									$homoPlus = $homo_1;
								}
								else{
									$homoPlus = $homo_2;
								}
								if($classe eq $homoPlus){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $homoPlus) and ($classe ne "NN")){
									print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($classe eq "$femelleb1$maleb1") or ($classe eq "$maleb1$femelleb1"))){
									print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "\n";
						}
						#obs: NN x AC real: dm null allele: A0 x AC
						elsif(($femelle eq "NN") and ($maleb1 ne $maleb2)){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<a0xab>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								
								my $homo_1 = $clasHOMO[0];
								my $homo_2 = $clasHOMO[1];
								
								my ($homoPlus, $homoMoins);
								
								if(($indvb1 eq $indvb2) and ($classe ne "NN") and ($HOMO_count{$homo_1} > $HOMO_count{$homo_2})){
									$homoPlus = $homo_1;
								}
								else{
									$homoPlus = $homo_2;
								}
								
								if($classe eq $male){
									print OUT "ab\t";
								}
								elsif($classe eq $homoPlus){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $homoPlus) and ($classe ne "NN")){
									print OUT "b0\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "\n";
						}
						#obs: AC x NN real: null allele: AC x A0
						elsif(($femelleb1 ne $femelleb2) and ($male eq "NN")){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxa0>\t";
							for(my $i = 6;$i<@tab;$i++){
								$tab[$i] =~ m/(\w)(\w)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = $tab[$i];
								
								my $homo_1 = $clasHOMO[0];
								my $homo_2 = $clasHOMO[1];
								
								my ($homoPlus, $homoMoins);
								
								if(($indvb1 eq $indvb2) and ($classe ne "NN") and ($HOMO_count{$homo_1} > $HOMO_count{$homo_2})){
									$homoPlus = $homo_1;
								}
								else{
									$homoPlus = $homo_2;
								}
								
								if($classe eq $femelle){
									print OUT "ab\t";
								}
								elsif($classe eq $homoPlus){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $homoPlus) and ($classe ne "NN")){
									print OUT "b0\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
							print OUT "\n";
						}
					}
					else{
						print SUP "prblHomoHetero	$line\n";
					}
				}
				else{
					print SUP "hetero>1	$line\n";
				}
			}
			else{
				print SUP "cls>3	$line\n";
			}
		}
		# case : same observed genotype for femele and male 
		#obs: AA x AA 
		elsif(($male eq $femelle) and ($femelleb1 eq $femelleb2) and ($male ne "NN") and ($femelle ne "NN")){
				# case : aa x ab ou ab x aa
				#obs: AA x AA ,real: err: AA x AC ou AC x AA => double
			if ($nbr_classes == 2){
				#print OUT "$line\n";
				print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
				for(my $i = 6;$i<@tab;$i++){
					$tab[$i] =~ m/(\w)(\w)/;
					$indvb1 = $1;
					$indvb2 = $2;
					$classe = $tab[$i];
					if(($classe eq $femelle) and ($classe ne "NN")){
						print OUT "aa\t";
					}
					elsif($indvb1 ne $indvb2){
						print OUT "ab\t";
					}
					elsif($classe eq "NN"){
						print OUT "--\t";
					}
					else{
						print OUT "NO\t";
					}
				}
				print OUT "\n";
				print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxaa>\t";
				for(my $i = 6;$i<@tab;$i++){
					$tab[$i] =~ m/(\w)(\w)/;
					$indvb1 = $1;
					$indvb2 = $2;
					$classe = $tab[$i];
					if(($classe eq $femelle) and ($classe ne "NN")){
						print OUT "aa\t";
					}
					elsif($indvb1 ne $indvb2){
						print OUT "ab\t";
					}
					elsif($classe eq "NN"){
						print OUT "--\t";
					}
					else{
						print OUT "NO\t";
					}
				}
					print OUT "\n";
			}
			elsif($nbr_classes == 3){
				if($nbr_hetero == 1){
					if ($countHETERO <= ($countHOMO/2)){
					# case : a0xab ou ab x a0 (
					#obs: AA x AA ,real: err null allele: A0 x AC ou AC x A0 => double
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<a0xab>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $femelle){
								print OUT "a-\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
								print OUT "b0\t";
							}
							elsif($indvb1 ne $indvb2){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxa0>\t";
						for(my $i = 6;$i<@tab;$i++){
							$tab[$i] =~ m/(\w)(\w)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = $tab[$i];
							if($classe eq $femelle){
								print OUT "a-\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
								print OUT "b0\t";
							}
							elsif($indvb1 ne $indvb2){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "\n";
					}
					else{
						print SUP "AAxAAprblHomoHetero	$line\n";
					}
				}
				else{
					print SUP "AAxAAHetro>1	$line\n";
				}
			}
			else{
				print SUP "AAxAAPlus3clas	$line\n";
			}
		}
		else{
			print SUP "2NN	$line\n";
		}
	}
	else{
		print SUP "cls=1	$line\n";
	}
	}
	
$nbr_classes = 0;
%CLASSES = ();
%HOMO = ();
%HETERO = ();
@clasHOMO = ();
%count = ();
%HOMO_count  = ();

$countHOMO = 0;
$countHETERO = 0;

$nbr_hetero = 0;
$nbr_homo = 0;

}

close F;
close OUT;


=head1 INCOMPATIBILITIES

Fully compatble with any perl version

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 2

=item Hajar CHOUIKI (INRA), hajar.chouiki-at-supagro.inra.fr

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


