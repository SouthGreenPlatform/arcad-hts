#!/bin/env perl

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

=pod

=head1 NAME

marker_to_ab.pl - transforms 01 genotype coding (phased vcf file) into nucleotidic coding (csv file) and assemble phased SNPs

=head1 DESCRIPTION

Transforms 01 genotype coding into nucleotidic coding and assemble phased SNPs

If a genotypic class is present only once, transform it into missing data

=head1 SYNOPSIS / USAGE

marker_to_ab.pl	[--help] [--man] [--version]
[-i Phased-VCF-file] \
[-o output-file] \


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

Phased VCF file


=item B<[-o]> ([output_file]): 

output file


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
if ($version) {print "1.0.0\n"; exit;}

my $out = $output.".csv";
my $sup = $output."_removed.csv";
my $info = $output."_info.csv";

open (F, $input)or die "cannot open the file $input\n";
open(OUT, ">$out") or die"cannot open the file $out\n";
open(SUP, ">$sup") or die"cannot open the file $sup\n";
open(INFO, ">$info") or die"cannot open the file $info\n";

my (@tab, @tab_base_alt, %HOMO_count, @not, @NotOk, @clas1, @clas2, @clasHETERO, @clasHOMO, @heterob1, @heterob2, @homob1, @homob2, @baseHetero, @baseHetero1, %hachage, $etat_marker, $baseCom, $alt_base1, $alt_base2, $alt_base3, $contig, $pos, $base_ref, $base_alt, $male, $femelle, $maleb1, $maleb2, $femelleb1, $femelleb2, $indvb1, $indvb2 );
my (%CLASSES, %HOMO, %HETERO, $indiv_marker, $group_marker, $marker, $countHOMO, $countHETERO, %nbr_HETERO, %nbr_HOMO, $classe, $nbr_classes, $total_indv, $nbr_hetero, $nbr_homo );
my ($_1,  $_2, $_3, $_4, $_5, $_6, $_7, $_8, $_9, $_10, $_11, $_12, $_13, $_14, $_15, $_16, $_17, $_18, $_19, $_20, $_21, $_22, $_23, $_24, $_25, $_26, $_27, $_28, $_29, $_30, $_31, $_32, $_33, $_34, $_35, $_36, $_37, $_38, $_39, $_40, $_41, $_42, $_43, $_44, $_45, $_46, $_47, $_48, $_49, $_50, $_51, $_52, $_53, $_54, $_55, $_56, $_57, $_58, $_59);
$_1= $_2= $_3= $_4= $_5= $_6= $_7= $_8= $_9= $_10= $_11= $_12= $_13= $_14= $_15= $_16= $_17= $_18= $_19= $_20= $_21= $_22= $_23= $_24= $_25= $_26= $_27= $_28= $_29= $_30= $_31= $_32= $_33= $_34= $_35= $_36= $_37= $_38= $_39= $_40= $_41= $_42= $_43= $_44= $_45= $_46= $_47= $_48= $_49= $_50= $_51= $_52= $_53= $_54= $_55= $_56= $_57= $_58= $_59 = 0;

$total_indv = 0;
$nbr_classes = 0;
$nbr_homo = 0;
$nbr_hetero = 0;
$countHOMO = 0;
$countHETERO = 0;

$indiv_marker = 0;
$group_marker = 0;

@clas1 = ();
@clas2 = ();
@clasHOMO = ();
@clasHETERO = ();
@heterob1 = ();
@heterob2 = ();
@homob1 = ();
@homob2 = ();
@baseHetero = ();
@baseHetero1 = ();
my %count = ( );
%hachage = ( );
%HOMO_count = ();
# Header :
print OUT "Chromosome	POS	Marker	Segregation	";

while(my $line =<F>){
	chomp $line;
	@tab = split (/\t/, $line);
	#la suite du header : les indvidus
	if($line =~ m/^Chromosome/){
		@tab = split (/\t/, $line);
		for(my $i = 4;$i<@tab;$i++){
			print OUT "$tab[$i]\t";
		}
		print OUT "\n";
	}
	else{
		$contig = $tab[0];
		$pos = $tab[1];
		$femelle = $tab[2];
		$male = $tab[3];
		$male =~ m/(\w+)\|(\w+)/;
		$maleb1 = $1;
		$maleb2 = $2;
		$femelle =~ m/(\w+)\|(\w+)/;
		$femelleb1 = $1;
		$femelleb2 = $2;
		
		# marqueur individuel : une seule base
		if( (length($femelleb1) == 1) and (length($maleb1) == 1)){
			my $femelleA = $tab[2];
			$femelleA =~ m/(\w+)\|(\w+)/;
			my $GNfemelleA = "$1$2";
			my @GNfemelleA = split ('', $GNfemelleA);
			my @chfemelleA = sort({$a cmp $b} @GNfemelleA );
			my $GNTfemelleA=join('', @chfemelleA);
			$femelle = "$GNTfemelleA";
			$femelle =~ m/(\w)(\w)/;
			$femelleb1 = $1;
			$femelleb2 = $2;
			
			my $maleA = $tab[3];
			$maleA =~ m/(\w+)\|(\w+)/;
			my $GNmaleA = "$1$2";
			my @GNmaleA = split ('', $GNmaleA);
			my @chmaleA = sort({$a cmp $b} @GNmaleA );
			my $GNTmaleA=join('', @chmaleA);
			$male = "$GNTmaleA";
			$male =~ m/(\w)(\w)/;
			$maleb1 = $1;
			$maleb2 = $2;
			$indiv_marker++;
			$etat_marker = "indiv_marker";
			
		}
		# marqueur regroupé : plusieurs bases
		else{
			$femelle = $tab[2];
			$male = $tab[3];
			$male =~ m/(\w+)\|(\w+)/;
			$maleb1 = $1;
			$maleb2 = $2;
			$male = "$maleb1$maleb2";
			$femelle =~ m/(\w+)\|(\w+)/;
			$femelleb1 = $1;
			$femelleb2 = $2;
			$femelle = "$femelleb1$femelleb2";
			$group_marker++;
			$etat_marker = "group_marker";
		}
	#calcul du nombre de classes génotypique et de l'effectif de chaque classe
	for(my $i = 4;$i<@tab;$i++){
		if($tab[$i] =~ m/(\w+)\|(\w+)/){
			$classe = "$1$2";
			if((length($1) == 1) and (length($2) == 1)){
				my $GN = "$1$2";
				my @GN = split ('', $GN);
				my @ch = sort({$a cmp $b} @GN );
				my $GNT=join('', @ch);
				$classe = "$GNT";
			}
			if($classe ne "NN"){
				if((!exists($CLASSES{"$1$2"})) and (!exists($CLASSES{"$2$1"}))){
					push(@clas1, $1);
					push(@clas2, $2);
					$CLASSES{$classe} = 1;
					$nbr_classes++;
					if($1 eq $2){
						# Mettre les classes homo dans le tableau @clasHOMO et chaque alléle dans un tableau : 
						## alléle 1 dans @homob1 et alléle 2 dans @homob2 
						push(@clasHOMO, $classe);
						push(@homob1, $1);
						push(@homob2, $2);
						$HOMO{$classe} = 1;
						$nbr_homo++;
						$HOMO_count{$classe}++;
					}
					else{
						push(@clasHETERO, $classe);
						push(@heterob1, $1);
						push(@heterob2, $2);
						
						# Cette condition va nous servir dans les cas : 36 et 39 pour récupérer les alléles différents des parents
						if((($1 ne $maleb1) and ($1 ne $maleb2)) or (($1 ne $femelleb1) and ($1 ne $femelleb2))){
							push (@baseHetero1, $1);
						}
						elsif((($2 ne $maleb1) and ($2 ne $maleb2)) or (($2 ne $femelleb1) and ($2 ne $femelleb2))){
							push(@baseHetero1, $2);
						}
						# éliminer les doublons dans @baseHetero1
						%hachage = map { $_, 1 } @baseHetero1;
						@baseHetero = keys %hachage;
						
						$HETERO{$classe} = 1;
						$nbr_hetero++;
						$HETERO{$classe}++;
					}
				}
				if((exists($CLASSES{"$1$2"})) or (exists($CLASSES{"$2$1"}))){
					$count{$classe}++;
					if($1 eq $2){
						$HOMO{$classe} = 1;
						$countHOMO++;
						$HOMO_count{$classe}++;
					}
					else{
						$HETERO{$classe} = 1;
						$countHETERO++;
						$HETERO{$classe}++;
					}
				}
			}
		}
	}

	#recuperer que les SNPs qui ont plus d'une classe genotypique dans la descendance
	if($nbr_classes > 1){
		#1# ab x ab
		if((($femelle eq $male) or ($femelle eq "$maleb1$maleb2") or ($femelle eq "$maleb2$maleb1")) and ($femelleb1 ne $femelleb2)){
			if($nbr_classes == 3){
				#print OUT "$line\n";
				print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
				for(my $i = 4;$i<@tab;$i++){
					$tab[$i] =~ m/(\w+)\|(\w+)/;
					$indvb1 = $1;
					$indvb2 = $2;
					$classe = "$indvb1$indvb2";
					if((length($1) == 1) and (length($2) == 1)){
						my $chaine = "$indvb1$indvb2";
						my $GN = "$1$2";
						my @GN = split ('', $GN);
						my @ch = sort({$a cmp $b} @GN );
						my $GNT=join('', @ch);
						$classe = "$GNT";
					}
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
				print OUT "	1	$etat_marker	\n";
				$_1++;
			}
			#on ne doit pas tomber sur ce cas => mettre dans un autre fichier
			elsif($nbr_classes == 2){
				print SUP "2clas:abxab	$line\n";
			}
			else{
				print SUP "#3;2 clas:abxab	$line\n";
			}
		}
		# aa x ab ou ab x aa
		elsif(($femelle ne $male) or (($femelle ne "$maleb1$maleb2") and ($femelle ne "$maleb2$maleb1"))){
			if ($nbr_classes == 2){
				#aa x ab
				#2# obs = reel : AA x AC ou AA x CA 
				if(($femelle ne "NN") and ($femelleb1 eq $femelleb2) and ($maleb1 ne $maleb2) and (($femelleb1 eq $maleb1) or ($femelleb1 eq $maleb2) or ($femelleb2 eq $maleb1) or ($femelleb2 eq $maleb2))){
					##print OUT "$line\n";
					print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
					for(my $i = 4;$i<@tab;$i++){
						$tab[$i] =~ m/(\w+)\|(\w+)/;
						$indvb1 = $1;
						$indvb2 = $2;
						$classe = "$indvb1$indvb2";
						
						if((length($1) == 1) and (length($2) == 1)){
							my $chaine = "$indvb1$indvb2";
							my $GN = "$1$2";
							my @GN = split ('', $GN);
							my @ch = sort({$a cmp $b} @GN );
							my $GNT=join('', @ch);
							$classe = "$GNT";
						}
						
						if((($classe eq $femelle) and ($indvb1 eq $indvb2)) and ($classe ne "NN")){
							print OUT "aa\t";
						}
						elsif(($indvb1 ne $indvb2) and (($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1") or ($classe eq $male))){
							print OUT "ab\t";
						}
						elsif($classe eq "NN"){
							print OUT "--\t";
						}
						# Juste pour vérifier s'il y a un génotype qui ne va pas, on le note "NO"
						else{
							print OUT "NO\t";
						}
					}
						print OUT "	2	$etat_marker	\n";
						$_2++;
				}
				#3# obs = reel : AA x CG
				elsif(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and ($femelleb1 ne $maleb1) and ($femelleb1 ne $maleb2) and ($maleb1 ne $maleb2)){
					if($nbr_homo == 0){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($indvb1 ne $indvb2) and (($classe eq "$femelleb1$maleb1") or ($classe eq "$maleb1$femelleb1"))){
								print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($classe eq "$femelleb1$maleb2") or ($classe eq "$maleb2$femelleb1"))){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	3	$etat_marker	\n";
						$_3++;
					}
					else{
						print SUP "2clas:obs=reel:AAxCG_homo#0	$line\n";
					}
				}
				#4# obs : AA x CC ou (CC x AA) ,reel : allele nul: A0 x CC ou (allele nul: CC x A0 ou err: AC x AA ou AC x AA)) => doubler
				elsif(($femelle ne "NN") and ($male ne "NN") and ($femelleb1 eq $femelleb2) and ($maleb1 eq $maleb2)){
					if($nbr_homo == 1){
						#print SUP "S:AAXCC_aaxab_err_A0xCC	$line\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($classe ne "NN") and ($indvb1 eq $indvb2)){
								#print OUT "aa\t";
							}
							elsif($indvb1 ne $indvb2){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						#print OUT "\n";
						##print OUT "$line\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
						if((length($1) == 1) and (length($2) == 1)){
							my $chaine = "$indvb1$indvb2";
							my $GN = "$1$2";
							my @GN = split ('', $GN);
							my @ch = sort({$a cmp $b} @GN );
							my $GNT=join('', @ch);
							$classe = "$GNT";
						}
							
							if(($classe ne "NN") and ($indvb1 eq $indvb2)){
								#print OUT "aa\t";
							}
							elsif($indvb1 ne $indvb2){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
							print SUP "err:4	$line\n";
							$_4++;
					}
					#5# obs: AA X CC ou CC x AA, err, reel: AA x CG ou CG x AA => doubler
					elsif($nbr_homo == 0){
						#print SUP "S:AAXCC_aaxab_err_AAxCG	$line\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($indvb1 ne $indvb2) and (($classe eq "$maleb1$femelleb1") or ($classe eq "$femelleb1$maleb1"))){
								#print OUT "ab\t";
							}
							elsif(($indvb1 ne $indvb2) and ($classe ne "$maleb1$femelleb1") and ($classe ne "$femelleb1$maleb1")){
								#print OUT "aa\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						#print OUT "\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($indvb1 ne $indvb2) and (($classe eq "$maleb1$femelleb1") or ($classe eq "$femelleb1$maleb1"))){
								#print OUT "ab\t";
							}
							elsif(($indvb1 ne $indvb2) and ($classe ne "$maleb1$femelleb1") and ($classe ne "$femelleb1$maleb1")){
								#print OUT "aa\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "err:5	$line\n";
						$_5++;
					}
					else{
						print SUP "2clas:AAxCC:HOMO#0,1	$line\n";
					}
				}
				# obs: AA x NN
				elsif(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and ($male eq "NN")){
					#6# obs: AA x NN, reel :(dm: AA x AC ou dm allele nul: CC x A0 ou A0 x CC ou err+dm: AC x AA) => doubler
					if($nbr_hetero == 1){
						#print OUT "$line\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
								#print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq $femelleb1) or ($indvb2 eq $femelleb1))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
							#print OUT "\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxaa>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
								
								if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
									#print OUT "aa\t";
								}
								elsif(($indvb1 ne $indvb2) and (($indvb1 eq $femelleb1) or ($indvb2 eq $femelleb1))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
						print SUP "err:6	$line\n";
						$_6++;
					}
					#7# AA x NN ,reel: err allele nul:AC x 00
					elsif($nbr_hetero == 0){
						#print OUT "$line\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
								#print OUT "aa\t";
							}
							elsif(($indvb1 eq $indvb2) and (($classe ne "$femelleb2$femelleb1") and ($classe ne "$femelleb1$femelleb2")) and ($classe ne "NN")){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "err:7	$line\n";
						$_7++;
					}
					#8# obs: AA x NN , dm reel: AA x CG ou err+dm reel:CG x AA => doubler
					elsif($nbr_homo == 0){
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
								my $all1_1 = $clas1[0];
								my $all1_2 = $clas2[0];
								my $all2_1 = $clas1[1];
								my $all2_2 = $clas2[1];
								my ($bas1, $bas2);
								if($all1_1 ne $femelleb1){
									$bas1 = $all1_1;
								}
								else{
									$bas1 = $all1_2;
								}
								if($all2_1 ne $femelleb1){
									$bas2 = $all2_2;
								}
								else{
									$bas2 = $all2_1;
								}
								
							if(($indvb1 ne $indvb2) and (($indvb1 eq "$femelleb1") or ($indvb2 eq "$femelleb1")) and (($indvb1 eq $bas1) or ($indvb2 eq $bas1))){
								#print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq "$femelleb1") or ($indvb2 eq "$femelleb1")) and (($indvb1 eq $bas2) or ($indvb2 eq $bas2))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						#print OUT "\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
								my $all1_1 = $clas1[0];
								my $all1_2 = $clas2[0];
								my $all2_1 = $clas1[1];
								my $all2_2 = $clas2[1];
								my ($bas1, $bas2);
								if($all1_1 ne $femelleb1){
									$bas1 = $all1_1;
								}
								else{
									$bas1 = $all1_2;
								}
								if($all2_1 ne $femelleb1){
									$bas2 = $all2_2;
								}
								else{
									$bas2 = $all2_1;
								}
								
							if(($indvb1 ne $indvb2) and (($indvb1 eq "$femelleb1") or ($indvb2 eq "$femelleb1")) and (($indvb1 eq $bas1) or ($indvb2 eq $bas1))){
								#print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq "$femelleb1") or ($indvb2 eq "$femelleb1")) and (($indvb1 eq $bas2) or ($indvb2 eq $bas2))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "err:8	$line\n";
						$_8++;
					}
					else{
						print SUP "2clas:AAxNN:HETERO>1;HOMO#0	$line\n";
					}
				}
				# obs: NN x AC
				elsif(($femelle eq "NN") and ($maleb1 ne $maleb2)){
					#9# obs: NN x AC reel:dm:  AA x AC
					if ($nbr_hetero == 1){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
						if((length($1) == 1) and (length($2) == 1)){
							my $chaine = "$indvb1$indvb2";
							my $GN = "$1$2";
							my @GN = split ('', $GN);
							my @ch = sort({$a cmp $b} @GN );
							my $GNT=join('', @ch);
							$classe = "$GNT";
						}
							
							if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
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
						print OUT "	9	$etat_marker	\n";
						$_9++;
					}
					#10# obs: NN X AC ,reel:allele nul :00 x AC
					elsif($nbr_hetero == 0){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
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
						print OUT "	10	$etat_marker	\n";
						$_10++;
					}
					#11# obs: NN  x CG, dm, reel: AA x CG
					elsif($nbr_homo == 0){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($indvb1 ne $indvb2) and (($indvb1 eq "$maleb1") and ($indvb2 ne "$maleb2")) or (($indvb2 eq "$maleb1") and ($indvb1 ne "$maleb2"))){
								print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq "$maleb2") and ($indvb2 ne "$maleb1")) or (($indvb2 eq "$maleb2") and ($indvb1 ne "$maleb1"))){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	11	$etat_marker	\n";
						$_11++;
					}
					else{
						print SUP "2clas:NNxAC:aaxab:hetero#0/#1;homo#0	$line\n";
					}
				}
				# ab x aa ou aa x ab
				#12# obs = reel : AC x AA 
				elsif(($male ne "NN") and ($femelle ne "NN") and ($maleb1 eq $maleb2) and ($femelleb1 ne $femelleb2) and (($maleb1 eq $femelleb1) or ($maleb1 eq $femelleb2))){
					#print OUT "$line\n";
					print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
					for(my $i = 4;$i<@tab;$i++){
						$tab[$i] =~ m/(\w+)\|(\w+)/;
						$indvb1 = $1;
						$indvb2 = $2;
						$classe = "$indvb1$indvb2";
						
						if((length($1) == 1) and (length($2) == 1)){
							my $chaine = "$indvb1$indvb2";
							my $GN = "$1$2";
							my @GN = split ('', $GN);
							my @ch = sort({$a cmp $b} @GN );
							my $GNT=join('', @ch);
							$classe = "$GNT";
						}
						
						if(($classe eq $male) and ($classe ne "NN" )){
							print OUT "aa\t";
						}
						elsif(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
							print OUT "ab\t";
						}
						elsif($classe eq "NN"){
							print OUT "--\t";
						}
						else{
							print OUT "NO\t";
						}
					}
						print OUT "	12	$etat_marker	\n";
						$_12++;
				}
				#13# obs = reel : CG x AA
				elsif(($maleb1 eq $maleb2) and ($male ne "NN") and ($maleb1 ne $femelleb1) and ($maleb1 ne $femelleb2) and ($femelleb1 ne $femelleb2)){
					if($nbr_homo == 0){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($indvb1 ne $indvb2) and (($classe eq "$maleb1$femelleb1") or ($classe eq "$femelleb1$maleb1"))){
								print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($classe eq "$maleb1$femelleb2") or ($classe eq "$femelleb2$maleb1"))){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	13	$etat_marker	\n";
						$_13++;
					}
					else{
						print SUP "2clas:HOMO>0:abxaa:CGxAA	$line\n";
					}
				}
				# NN x AA 
				elsif(($maleb1 eq $maleb2) and ($femelle eq "NN") and ($male ne "NN")){
					#14# obs: NN x AA (reel dm: AC x AA ou dm:allele nul: CC x A0 ou A0 x CC ou err+dm: AA x AC) => doubler
					if($nbr_hetero == 1){
						#print OUT "$line\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
						if((length($1) == 1) and (length($2) == 1)){
							my $chaine = "$indvb1$indvb2";
							my $GN = "$1$2";
							my @GN = split ('', $GN);
							my @ch = sort({$a cmp $b} @GN );
							my $GNT=join('', @ch);
							$classe = "$GNT";
						}
							
							if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
								#print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq $maleb1 ) or ($indvb2 eq $maleb1))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						#print OUT "\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
						if((length($1) == 1) and (length($2) == 1)){
							my $chaine = "$indvb1$indvb2";
							my $GN = "$1$2";
							my @GN = split ('', $GN);
							my @ch = sort({$a cmp $b} @GN );
							my $GNT=join('', @ch);
							$classe = "$GNT";
						}
							
							if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
								#print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq $maleb1 ) or ($indvb2 eq $maleb1))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "	err:14	\n";
						$_14++;
					}
					#15# NN x AA ,reel:  err+dm ,err allele nul:  00 x AC
					elsif($nbr_hetero == 0){
						#print OUT "$line\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
								#print OUT "aa\t";
							}
							elsif(($indvb1 eq $indvb2) and (($classe ne $male) and ($classe ne "$maleb1$maleb2") and ($classe ne "$maleb2$maleb1")) and ($classe ne "NN")){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "	err:15	\n";
						$_15++;
					}
					#16# obs: NN x AA, dm reel: CG x AA ou err+dm, reel: AA x CG => doubler
					elsif($nbr_homo == 0){
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
								my $all1_1 = $clas1[0];
								my $all1_2 = $clas2[0];
								my $all2_1 = $clas1[1];
								my $all2_2 = $clas2[1];
								my ($bas1, $bas2);
								if($all1_1 ne $maleb1){
									$bas1 = $all1_1;
								}
								else{
									$bas1 = $all1_2;
								}
								if($all2_1 ne $maleb1){
									$bas2 = $all2_2;
								}
								else{
									$bas2 = $all2_1;
								}
							
							if(($indvb1 ne $indvb2) and (($indvb1 eq "$maleb1") or ($indvb2 eq "$maleb1")) and (($indvb1 eq $bas1) or ($indvb2 eq $bas1))){
								#print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq "$maleb1") or ($indvb2 eq "$maleb1")) and (($indvb1 eq $bas2) or ($indvb2 eq $bas2))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						#print OUT "\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<aaxab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
								my $all1_1 = $clas1[0];
								my $all1_2 = $clas2[0];
								my $all2_1 = $clas1[1];
								my $all2_2 = $clas2[1];
								my ($bas1, $bas2);
								if($all1_1 ne $maleb1){
									$bas1 = $all1_1;
								}
								else{
									$bas1 = $all1_2;
								}
								if($all2_1 ne $maleb1){
									$bas2 = $all2_2;
								}
								else{
									$bas2 = $all2_1;
								}
							
							if(($indvb1 ne $indvb2) and (($indvb1 eq "$maleb1") or ($indvb2 eq "$maleb1")) and (($indvb1 eq $bas1) or ($indvb2 eq $bas1))){
								#print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq "$maleb1") or ($indvb2 eq "$maleb1")) and (($indvb1 eq $bas2) or ($indvb2 eq $bas2))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "	err:16	\n";
						$_16++;
					}
				}
				# AC x NN 
				elsif(($male eq "NN") and ($femelleb1 ne $femelleb2)){
					#17# obs: AC x NN reel:dm: AC x AA
					if($nbr_hetero == 1){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
								print OUT "ab\t";
							}
							elsif(($indvb1 eq $indvb2) and (($indvb1 eq $femelleb1) or ($indvb1 eq $femelleb2)) and ($classe ne "NN")){
								print OUT "aa\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
							print OUT "	17	$etat_marker	\n";
							$_17++;
					}
					#18# obs: AC X NN ,reel:allele nul: AC x 00
					elsif($nbr_hetero == 0){
						#print OUT "$line\n";
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
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
						print OUT "	18	$etat_marker	\n";
						$_18++;
					}
					#19# obs: CG x NN, dm, reel: CG x AA
					elsif($nbr_homo == 0){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxaa>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($indvb1 ne $indvb2) and (($indvb1 eq "$femelleb1") and ($indvb2 ne "$femelleb2")) or (($indvb2 eq "$femelleb1") and ($indvb1 ne "$femelleb2"))){
								print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq "$femelleb2") and ($indvb2 ne "$femelleb1")) or (($indvb2 eq "$femelleb2") and ($indvb1 ne "$femelleb1"))){
								print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	19	$etat_marker	\n";
						$_19++;
					}
					else{
						print SUP "2clas:abxaa:ACxNN:HETERO#0,1;HOMO#0	$line\n";
					}
				}
			}
			# abxab 
			elsif($nbr_classes == 3){
				if($nbr_hetero == 1){
					# abxab
					if($countHETERO >= 1.5 * ($countHOMO/2)){
						#20# obs: NN x AC reel:dm:  AC x AC
						if(($femelle eq "NN") and ($maleb1 ne $maleb2)){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
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
								print OUT "	20	$etat_marker	\n";
								$_20++;
						}
						#21# obs: AC x NN reel:dm: AC x AC
						elsif(($male eq "NN") and ($femelleb1 ne $femelleb2)){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
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
								print OUT "	21	$etat_marker	\n";
								$_21++;
						}
						#22# obs: AC x AA err, reel: AC x AC
						elsif(($femelleb1 ne $femelleb2) and ($maleb1 eq $maleb2) and ($male ne "NN") and (($maleb1 eq $femelleb1) or ($maleb1 eq $femelleb2))){
							#print OUT "$line\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
									#print OUT "ab\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe eq $male)){
									#print OUT "aa\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne "$male") and ($classe ne "NN")){
									#print OUT "bb\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
								print SUP "err:22	$line\n";
								$_22++;
						}
						#23# obs: AA x AC reel: err, AC x AC
						elsif(($femelleb1 eq $femelleb2) and ($maleb1 ne $maleb2) and ($femelle ne "NN") and (($femelleb1 eq $maleb1) or ($femelleb1 eq $maleb2))){
							#print OUT "$line\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
									#print OUT "ab\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe eq $femelle)){
									#print OUT "aa\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
									#print OUT "bb\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
								print SUP "err:23	$line\n";
								$_23++;
						}
						#24# obs: NN x AA, dm+err reel: AC x AC
						elsif(($femelle eq "NN") and ($maleb1 eq $maleb2) and ($male ne "NN")){
							#print OUT "$line\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $male) and ($classe ne "NN")){
									#print OUT "aa\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $male) and ($classe ne "NN")){
									#print OUT "bb\t";
								}
								elsif(($indvb1 ne $indvb2) and (($maleb1 ne $indvb1) or ($maleb1 ne $indvb2))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
								print SUP "err:24	$line\n";
								$_24++;
						}
						#25# obs: AA x NN, dm+err reel: AC x AC
						elsif(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and ($male eq "NN")){
							#print OUT "$line\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $femelle) and ($classe ne "NN")){
									#print OUT "aa\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne "NN") and ($classe ne $femelle)){
									#print OUT "bb\t";
								}
								elsif(($indvb1 ne $indvb2) and (($femelleb1 ne $indvb1) or ($femelleb1 ne $indvb2))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
								print SUP "err:25	$line\n";
								$_25++;
						}
						else{
							print SUP "3clas:abxab:AAxCC	$line\n";
						}
					}
					# a0xab or abxa0
					elsif ($countHETERO <= ($countHOMO/2)){
						#26# obs : AA X AC ou (AA x NN), reel: allele nul: A0 x AC ou (dm allele nul: A0 x AC ou err+dm allele nul: AC x A0) => doubler
						if(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and (($male eq "NN") or (($maleb1 ne $maleb2) and (($femelleb1 eq $maleb1) or ($femelleb1 eq $maleb2))))){
							#print OUT "$line\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<a0xab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
									#print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and (($classe ne $femelle) or ($classe ne "$femelleb1$femelleb2") or ($classe ne "$femelleb2$femelleb1")) and ($classe ne "NN")){
									#print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($femelleb1 eq $indvb1) or ($femelleb1 eq $indvb2))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
							#print OUT "\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxa0>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
									#print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and (($classe ne $femelle) or ($classe ne "$femelleb1$femelleb2") or ($classe ne "$femelleb2$femelleb1")) and ($classe ne "NN")){
									#print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($femelleb1 eq $indvb1) or ($femelleb1 eq $indvb2))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
								print SUP "err:26	$line\n";
								$_26++;
						}
						#27# obs : AC X AA ou (NN x AA), reel: allele nul: AC x A0 ou (dm allele nul: AC x A0 ou err+dm allele nul: A0 x AC) => doubler
						elsif(($maleb1 eq $maleb2) and ($male ne "NN") and (($femelle eq "NN") or (($femelleb1 ne $femelleb2) and (($maleb1 eq $femelleb1) or ($maleb1 eq $femelleb2))))){
							#print OUT "$line\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxa0>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
									#print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and (($classe ne $male) or ($classe ne "$maleb1$maleb2") or ($classe ne "$maleb2$maleb1")) and ($classe ne "NN")){
									#print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($maleb1 eq $indvb1) or ($maleb1 eq $indvb2))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
							#print OUT "\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<a0xab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
								if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
									#print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and (($classe ne $male) or ($classe ne "$maleb1$maleb2") or ($classe ne "$maleb2$maleb1")) and ($classe ne "NN")){
									#print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($maleb1 eq $indvb1) or ($maleb1 eq $indvb2))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
								print SUP "err:27	$line\n";
								$_27++;
						}
						#28# obs: AA x CC ou CC x AA, reel: err: A0 x AC ou AC x A0 => doubler
						elsif(($femelleb1 eq $femelleb2) and ($maleb1 eq $maleb2) and ($male ne "NN") and ($femelle ne "NN")){
							#print OUT "$line\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<a0xab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								
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
									#print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne "NN") and ($classe ne $homoPlus)){
									#print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($classe eq "$femelleb1$maleb1") or ($classe eq "$maleb1$femelleb1"))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
							#print OUT "\n";
							#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxa0>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								##!!! homo le plus frequent 
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
									#print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne "NN") and ($classe ne $homoPlus)){
									#print OUT "b0\t";
								}
								elsif(($indvb1 ne $indvb2) and (($classe eq "$femelleb1$maleb1") or ($classe eq "$maleb1$femelleb1"))){
									#print OUT "ab\t";
								}
								elsif($classe eq "NN"){
									#print OUT "--\t";
								}
								else {
									#print OUT "NO\t";
								}
							}
								print SUP "err:28	$line\n";
								$_28++;
						}
						#29# obs: NN x AC reel: dm allele nul: A0 x AC
						elsif(($femelle eq "NN") and ($maleb1 ne $maleb2)){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<a0xab>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								my $homo_1 = $clasHOMO[0];
								my $homo_2 = $clasHOMO[1];
								
								my ($homoPlus, $homoMoins);
								
								if(($indvb1 eq $indvb2) and ($classe ne "NN") and ($HOMO_count{$homo_1} > $HOMO_count{$homo_2})){
									$homoPlus = $homo_1;
								}
								else{
									$homoPlus = $homo_2;
								}
								if(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
									print OUT "ab\t";
								}
								# homo le + frequent
								elsif($classe eq $homoPlus){
									print OUT "a-\t";
								}
								# homo le - frequent
								elsif(($indvb1 eq $indvb2) and ($classe ne $homoPlus) and ($classe ne "NN") ){
									print OUT "b0\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
								print OUT "	29	$etat_marker	\n";
								$_29++;
						}
						#30# obs: AC x NN reel: allele nul: AC x A0
						elsif(($male eq "NN") and ($femelleb1 ne $femelleb2)){
							#print OUT "$line\n";
							print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxa0>\t";
							for(my $i = 4;$i<@tab;$i++){
								$tab[$i] =~ m/(\w+)\|(\w+)/;
								$indvb1 = $1;
								$indvb2 = $2;
								$classe = "$indvb1$indvb2";
								
								if((length($1) == 1) and (length($2) == 1)){
									my $chaine = "$indvb1$indvb2";
									my $GN = "$1$2";
									my @GN = split ('', $GN);
									my @ch = sort({$a cmp $b} @GN );
									my $GNT=join('', @ch);
									$classe = "$GNT";
								}
								my $homo_1 = $clasHOMO[0];
								my $homo_2 = $clasHOMO[1];
								
								my ($homoPlus, $homoMoins);
								
								if(($indvb1 eq $indvb2) and ($classe ne "NN") and ($HOMO_count{$homo_1} > $HOMO_count{$homo_2})){
									$homoPlus = $homo_1;
								}
								else{
									$homoPlus = $homo_2;
								}
								# homo le plus frequent
								if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
									print OUT "ab\t";
								}
								elsif($classe eq $homoPlus){
									print OUT "a-\t";
								}
								elsif(($indvb1 eq $indvb2) and ($classe ne $homoPlus) and ($classe ne "NN") ){
									print OUT "b0\t";
								}
								elsif($classe eq "NN"){
									print OUT "--\t";
								}
								else {
									print OUT "NO\t";
								}
							}
							print OUT "	30	$etat_marker	\n";
							$_30++;
						}
						else{
							print SUP "3clas:a0xabORabxa0=autre	$line\n";
						}
					}
					else{
						print SUP "3clas:prblHomoHetero	$line\n";
					}
				}
				else{
					print SUP "3clas:hetero>1	$line\n";
				}
			}
			elsif($nbr_classes == 4){
				# abxcd
				if(($femelleb1 ne $femelleb2) and ($maleb1 ne $maleb2) and ($femelleb1 ne $maleb1) and ($femelleb1 ne $maleb2) and ($femelleb2 ne $maleb1) and ($femelleb2 ne $maleb2)){
					#31# obs = reel : AC x GT 
					if($nbr_homo == 0){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($classe eq "$maleb1$femelleb1") or ($classe eq "$femelleb1$maleb1")){
								print OUT "ac\t";
							}
							elsif(($classe eq "$maleb1$femelleb2") or ($classe eq "$femelleb2$maleb1")){
								print OUT "bc\t";
							}
							elsif(($classe eq "$maleb2$femelleb1") or ($classe eq "$femelleb1$maleb2")){
								print OUT "ad\t";
							}
							elsif(($classe eq "$maleb2$femelleb2") or ($classe eq "$femelleb2$maleb2")){
								print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	31	$etat_marker	\n";
						$_31++;
					}
					#
					else{
						print SUP "3clas:abxcdreel_HOMO#0	$line\n";
					}
				}
				# obs: CG x AA 
				elsif(($femelleb1 ne $femelleb2) and ($maleb1 eq $maleb2) and ($male ne "NN") and ($maleb1 ne $femelleb1) and ($maleb1 ne $femelleb2)){
					#32# obs : CG x AA, all. nul: reel : CG x A0
					if($nbr_homo == 2){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($classe eq "$maleb1$femelleb1") or ($classe eq "$femelleb1$maleb1")){
								print OUT "ac\t";
							}
							elsif(($classe eq "$maleb1$femelleb2") or ($classe eq "$femelleb2$maleb1")){
								print OUT "bc\t";
							}
							elsif($classe eq "$femelleb1$femelleb1"){
								print OUT "ad\t";
							}
							elsif($classe eq "$femelleb2$femelleb2"){
								print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	32	$etat_marker	\n";
						$_32++;
					}
					#33#  obs : AC x GG, err: reel : AC x GT 
					elsif($nbr_homo == 0){
						print SUP "err:33	$line\n";
						$_33++;
					}
					else{
						print SUP "4clas:abxcd:AAxCG_HOMO != 0,2	$line\n";
					}
				}
				# obs: AA x CG 
				elsif(($maleb1 ne $maleb2) and ($femelleb1 eq $femelleb2) and ($femelle ne "NN") and ($femelleb1 ne $maleb1) and ($femelleb1 ne $maleb2)){
					#34# obs :AA x CG , all. nul: reel : A0 x CG 
					if($nbr_homo == 2){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($classe eq "$femelleb1$maleb1") or ($classe eq "$maleb1$femelleb1")){
								print OUT "ac\t";
							}
							elsif(($classe eq "$femelleb1$maleb2") or ($classe eq "$maleb2$femelleb1")){
								print OUT "ad\t";
							}
							elsif($classe eq "$maleb1$maleb1"){
								print OUT "bc\t";
							}
							elsif($classe eq "$maleb2$maleb2"){
								print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	34	$etat_marker	\n";
						$_34++;
					}
					#35# obs :AA x GT , err: reel : AC x GT 
					elsif($nbr_homo == 0){
						print SUP "err:35	$line\n";
						$_35++;
					}
					else{
						print SUP "4clas:abxcd:AAxGT_HOMO!=0,2	$line\n";
					}
				}
				# obs :NN x GT 
				elsif(($femelle eq "NN") and ($maleb1 ne $maleb2)){
					#36# obs :NN x GT , dm: reel : AC x GT 
					if($nbr_homo == 0){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							my ($bas1 ,$bas2);
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
								my $base1 = $baseHetero[0];
								my $base2 = $baseHetero[1];
								 
							if((($indvb1 eq $maleb1) and ($indvb2 eq $base1)) or (($indvb2 eq $maleb1) and ($indvb1 eq $base1))){
								print OUT "ac\t";
							}
							elsif((($indvb1 eq $maleb1) and ($indvb2 eq $base2)) or (($indvb2 eq $maleb1) and ($indvb1 eq $base2))){
								print OUT "bc\t";
							}
							elsif((($indvb1 eq $maleb2) and ($indvb2 eq $base1)) or (($indvb2 eq $maleb2) and ($indvb1 eq $base1))){
								print OUT "ad\t";
							}
							elsif((($indvb1 eq $maleb2) and ($indvb2 eq $base2)) or (($indvb2 eq $maleb2) and ($indvb1 eq $base2))){
								print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	36	$etat_marker	\n";
						$_36++;
					}
					# ab x ac
					#37# obs : NN x AG , dm: reel : AC x AG 
					elsif($nbr_homo == 1){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxac>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							my $all_com;
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							$all_com = $homob1[0];
							if(($indvb1 eq $indvb1) and ($classe ne "NN")){
								$all_com = $indvb1;
								print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb1) and ((($indvb1 eq $all_com) and (($indvb2 ne $maleb1) and ($indvb2 ne $maleb2))) or (($indvb2 eq $all_com) and ($indvb1 ne $maleb1) and ($indvb1 ne $maleb2)))){
								print OUT "ab\t";
							}
							elsif(($indvb1 ne $indvb1) and ($classe eq "$maleb1$maleb2") and ($classe eq "$maleb2$maleb1")){
								print OUT "ac\t";
							}
							elsif(($indvb1 ne $indvb1) and ($indvb1 ne $all_com) and ($indvb2 ne $all_com) and (($indvb1 eq $maleb2) or ($indvb2 eq $maleb2) or ($indvb1 eq $maleb1) or ($indvb2 eq $maleb1))){
								print OUT "bc\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	37	$etat_marker	\n";
						$_37++;
					}
					#38# obs : NN x CG, dm+all. nul reel: A0 x CG
					elsif($nbr_homo == 2){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $maleb1) and ($indvb2 ne $maleb2)) or (($indvb2 eq $maleb1) and ($indvb1 ne $maleb2)))){
								print OUT "ac\t";
							}
							elsif(($indvb1 ne $indvb2) and ((($indvb1 eq $maleb2) and ($indvb2 ne $maleb1)) or (($indvb2 eq $maleb2) and ($indvb1 ne $maleb1)))){
								print OUT "ad\t";
							}
							elsif($classe eq "$maleb1$maleb1"){
								print OUT "bc\t";
							}
							elsif($classe eq "$maleb2$maleb2"){
								print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	38	$etat_marker	\n";
						$_38++;
					}
				}
				# obs: AC x NN
				elsif(($male eq "NN") and ($femelleb1 ne $femelleb2)){
					#39# obs : AC x NN, dm: reel : AC x GT
					if($nbr_homo == 0){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
								my $base1 = $baseHetero[0];
								my $base2 = $baseHetero[1];
								
							if((($indvb1 eq $femelleb1) and ($indvb2 eq $base1)) or (($indvb2 eq $femelleb1) and ($indvb1 eq $base1))){
								print OUT "ac\t";
							}
							elsif((($indvb1 eq $femelleb1) and ($indvb2 eq $base2)) or (($indvb2 eq $femelleb1) and ($indvb1 eq $base2))){
								print OUT "ad\t";
							}
							elsif((($indvb1 eq $femelleb2) and ($indvb2 eq $base1)) or (($indvb2 eq $femelleb2) and ($indvb1 eq $base1))){
								print OUT "bc\t";
							}
							elsif((($indvb1 eq $femelleb2) and ($indvb2 eq $base2)) or (($indvb2 eq $femelleb2) and ($indvb1 eq $base2))){
								print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	39	$etat_marker	\n";
						$_39++;
					}
					# ab x ac
					#40# obs : AC x NN , dm: reel : AC x AG 
					elsif($nbr_homo == 1){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxac>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							my $all_com;
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							$all_com = $homob1[0];
							if(($indvb1 eq $indvb2) and ($classe ne "NN")){
								$all_com = $indvb1;
								print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and ((($indvb1 eq $all_com) and (($indvb2 ne $femelleb1) and ($indvb2 ne $femelleb2))) or (($indvb2 eq $all_com) and ($indvb1 ne $femelleb1) and ($indvb1 ne $femelleb2)))){
								print OUT "ac\t";
							}
							elsif(($indvb1 ne $indvb2) and ($classe eq "$femelleb1$femelleb2") and ($classe eq "$femelleb2$femelleb1")){
								print OUT "ab\t";
							}
							elsif(($indvb1 ne $indvb2) and ($indvb1 ne $all_com) and ($indvb2 ne $all_com) and (($indvb1 eq $femelleb2) or ($indvb2 eq $femelleb2) or ($indvb1 eq $femelleb1) or ($indvb2 eq $femelleb1))){
								print OUT "bc\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	40	$etat_marker	\n";
						$_40++;
					}
					#41# obs : CG x NN, dm+all. nul: reel : CG x A0
					elsif($nbr_homo == 2){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $femelleb1) and ($indvb2 ne $femelleb2)) or (($indvb2 eq $femelleb1) and ($indvb1 ne $femelleb2)))){
								print OUT "ac\t";
							}
							elsif(($indvb1 ne $indvb2) and ((($indvb1 eq $femelleb2) and ($indvb2 ne $femelleb1)) or (($indvb2 eq $femelleb2) and ($indvb1 ne $femelleb1)))){
								print OUT "bc\t";
							}
							elsif($classe eq "$femelleb1$femelleb1"){
								print OUT "ad\t";
							}
							elsif($classe eq "$femelleb2$femelleb2"){
								print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	41	$etat_marker	\n";
						$_41++;
					}
				}
				# obs: NN x GG 
				elsif(($femelle eq "NN") and ($maleb1 ne $maleb2)){
					#42# obs : NN x AA, err+dm: reel : AC x GT 
					if($nbr_homo == 0){
						print SUP "err:42	$line\n";
						$_42++;
					}
					#43# obs : NN x AA , err+dm: reel : AC x AG
					elsif($nbr_homo == 1){
						print SUP "err:43	$line\n";
						$_43++;
					}
					#44# obs : NN x AA, dm+all. nul: reel : CG x A0 / err+dm+all. nul: reel: A0 x CG => Doubler
					elsif($nbr_homo == 2){
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							 
							my ($homo1, $homo2);
							$homo1 = $clasHOMO[0];
							$homo2 = $clasHOMO[1];
							
							my ($homo1b1, $homo2b1);
							$homo1b1 = $homob1[0];
							$homo2b1 = $homob1[1];
							
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $maleb1) and ($indvb2 eq $homo1b1)) or (($indvb2 eq $maleb1) and ($indvb1 eq $homo1b1)))){
								#print OUT "ac\t";
							}
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $maleb1) and ($indvb2 eq $homo2b1)) or (($indvb2 eq $maleb1) and ($indvb1 eq $homo2b1)))){
								#print OUT "bc\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "$maleb1$maleb1") and ($classe ne "NN") and ($classe eq $homo1)){
								#print OUT "ad\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "$maleb2$maleb2") and ($classe ne "NN") and ($classe eq $homo2)){
								#print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						#print OUT "\n";
						$_44++;
						#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							 
							my ($homo1, $homo2);
							$homo1 = $clasHOMO[0];
							$homo2 = $clasHOMO[1];
							
							my ($homo1b1, $homo2b1);
							$homo1b1 = $homob1[0];
							$homo2b1 = $homob1[1];
							
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $maleb1) and ($indvb2 eq $homo1b1)) or (($indvb2 eq $maleb1) and ($indvb1 eq $homo1b1)))){
								#print OUT "ac\t";
							}
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $maleb1) and ($indvb2 eq $homo2b1)) or (($indvb2 eq $maleb1) and ($indvb1 eq $homo2b1)))){
								#print OUT "ad\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "$maleb1$maleb1") and ($classe ne "NN") and ($classe eq $homo1)){
								#print OUT "bc\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "$maleb2$maleb2") and ($classe ne "NN") and ($classe eq $homo2)){
								#print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "err:44	$line\n";
					}
				}
				# obs: AA x NN 
				elsif(($male eq "NN") and ($femelleb1 eq $femelleb2)){
					#45# obs : AA x NN, err+dm: reel : AC x GT
					if($nbr_homo == 0){
						print SUP "err:45	$line\n";
						$_45++;
					}
					# ab x ac
					#46# obs : AA x NN , err+dm: reel : AC x AG
					elsif($nbr_homo == 1){
						print SUP "err:46	$line\n";
						$_46++;
					}
					#47# obs : AA x NN, dm+all. nul: reel : A0 x CG / err+dm+all. nul: reel : CG x A0 => Doubler
					elsif($nbr_homo == 2){
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							my ($homo1, $homo2);
							$homo1 = $clasHOMO[0];
							$homo2 = $clasHOMO[1];
							
							my ($homo1b1, $homo2b1);
							$homo1b1 = $homob1[0];
							$homo2b1 = $homob1[1];
							
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $femelleb1) and ($indvb2 eq $homo1b1)) or (($indvb2 eq $femelleb1) and ($indvb1 eq $homo1b1)))){
								#print OUT "ac\t";
							}
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $femelleb1) and ($indvb2 eq $homo2b1)) or (($indvb2 eq $femelleb1) and ($indvb1 eq $homo2b1)))){
								#print OUT "ad\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "$femelleb1$femelleb1") and ($classe ne "NN") and ($classe eq $homo1)){
								#print OUT "bc\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "$femelleb2$femelleb2") and ($classe ne "NN") and ($classe eq $homo2)){
								#print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						#print OUT "\n";
						$_47++;
						#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxcd>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							my ($homo1, $homo2);
							$homo1 = $clasHOMO[0];
							$homo2 = $clasHOMO[1];
							
							my ($homo1b1, $homo2b1);
							$homo1b1 = $homob1[0];
							$homo2b1 = $homob1[1];
							
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $femelleb1) and ($indvb2 eq $homo1b1)) or (($indvb2 eq $femelleb1) and ($indvb1 eq $homo1b1)))){
								#print OUT "ac\t";
							}
							if(($indvb1 ne $indvb2) and ((($indvb1 eq $femelleb1) and ($indvb2 eq $homo2b1)) or (($indvb2 eq $femelleb1) and ($indvb1 eq $homo2b1)))){
								#print OUT "bc\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "$femelleb1$femelleb1") and ($classe ne "NN") and ($classe eq $homo1)){
								#print OUT "ad\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "$femelleb2$femelleb2") and ($classe ne "NN") and ($classe eq $homo2)){
								#print OUT "bd\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "err:47	$line\n";
					}
				}
				#48# obs : AA x CC , err+all. nul: reel : A0 x CG 
				elsif(($femelleb1 eq $femelleb2) and ($maleb1 eq $maleb2) and ($femelleb1 ne $maleb1) and ($male ne "NN") and ($femelle ne "NN")){
					print SUP "err:48	$line\n";
					$_48++;
				}
				# abxac
				elsif(($femelleb1 ne $femelleb2) and ($maleb1 ne $maleb2) and (($femelleb1 eq $maleb1) or ($femelleb1 eq $maleb2) or ($femelleb2 eq $maleb1) or ($femelleb2 eq $maleb2))){
					#49# obs = reel : AC x AG 
					
					if($nbr_homo == 1){
						print OUT "$contig\t$pos\t$contig"."_"."$pos\t<abxac>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
								
							}
							my $baseCom = $homob1[0];
							
							if(($classe eq $femelle) or ($classe eq "$femelleb1$femelleb2") or ($classe eq "$femelleb2$femelleb1")){
								print OUT "ab\t";
								
							}
							elsif(($classe eq $male) or ($classe eq "$maleb1$maleb2") or ($classe eq "$maleb2$maleb1")){
								print OUT "ac\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne "NN") and ((($indvb1 eq $maleb1) or ($indvb1 eq $maleb2)) and (($indvb1 eq $femelleb1) or ($indvb1 eq $femelleb2)))){
								$baseCom = $indvb1;
								print OUT "aa\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 ne $baseCom) and ($indvb2 ne $baseCom)) and ($classe ne "$femelleb1$femelleb2") and ($classe ne "$femelleb2$femelleb1") and ($classe ne "$maleb1$maleb2") and ($classe ne "$maleb2$maleb1") and (((($indvb1 eq $femelleb1) or ($indvb1 eq $femelleb2)) and (($indvb2 eq $maleb1) or ($indvb2 eq $maleb2))) or ((($indvb2 eq $femelleb1) or ($indvb2 eq $femelleb2)) and (($indvb1 eq $maleb1) or ($indvb1 eq $maleb2))))){
								print OUT "bc\t";
							}
							elsif($classe eq "NN"){
								print OUT "--\t";
							}
							else {
								print OUT "NO\t";
							}
						}
						print OUT "	49	$etat_marker	\n";
						$_49++;
					}
					else{
						print SUP "4clas:abxac_HOMO>1	$line\n";
					}
				}
				#50# AA X AG err, reel AC x AG 
				elsif(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and ($maleb1 ne $maleb2) and (($maleb1 eq $femelleb1) or ($maleb2 eq $femelleb1))){
					if($nbr_homo == 1){
					$_50++;
					print SUP "err:50	$line\n";
					}
					else{
						print SUP "4clas:AAxAG:#1HOMO	$line\n";
					}
				}
				#51# AC X AA err, reel AC x AG 
				elsif(($femelleb1 ne $femelleb2) and ($male ne "NN") and ($maleb1 eq $maleb2) and (($femelleb1 eq $maleb1) or ($femelleb2 eq $maleb1))){
					if($nbr_homo == 1){
					$_51++;
					print SUP "err:51	$line\n";
					}
					else{
						print SUP "4clas:ACxAA:#1HOMO	$line\n";
					}
				}
				#52# CC X AG err, reel AC x AG 
				elsif(($femelleb1 eq $femelleb2) and ($femelle ne "NN") and ($maleb1 ne $maleb2) and (($maleb1 ne $femelleb1) and ($maleb2 ne $femelleb1))){
					if($nbr_homo == 1){
					$_52++;
					print SUP ":err:52	$line\n";
					}
					else{
						print SUP "4clas:CCxAG:#1HOMO	$line\n";
					}
				}
				#53# AC X GG err, reel AC x AG 
				elsif(($femelleb1 ne $femelleb2) and ($male ne "NN") and ($maleb1 eq $maleb2) and (($femelleb1 ne $maleb1) and ($femelleb2 ne $maleb1))){
					if($nbr_homo == 1){
					$_53++;
					print SUP "err:53	$line\n";
					}
					else{
						print SUP "4clas:ACxGG:#1HOMO	$line\n";
					}
				}
				else{
					print SUP "4clas:autreCas	$line\n";
				}
			}
			else{
				print SUP "Plus4clas	$line\n";
			}
		}
		#obs: AA x AA 
		elsif(($male eq $femelle) and ($femelleb1 eq $femelleb2) and ($male ne "NN") and ($femelle ne "NN")){
				# aa x ab ou ab x aa
				#54# obs: AA x AA ,reel: err: AA x AC ou AC x AA => doubler 
			if ($nbr_classes == 2){
				#print SUP "S:AAXAAerrACXAA	$line\n"; 
				#print OUT "$line\n";
				#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<aaxab>\t";
				for(my $i = 4;$i<@tab;$i++){
					$tab[$i] =~ m/(\w+)\|(\w+)/;
					$indvb1 = $1;
					$indvb2 = $2;
					$classe = "$indvb1$indvb2";
					
					if((length($1) == 1) and (length($2) == 1)){
						my $chaine = "$indvb1$indvb2";
						my $GN = "$1$2";
						my @GN = split ('', $GN);
						my @ch = sort({$a cmp $b} @GN );
						my $GNT=join('', @ch);
						$classe = "$GNT";
					}
					
					if(($classe eq $femelle) and ($classe ne "NN")){
						#print OUT "aa\t";
					}
					elsif(($indvb1 ne $indvb2) and (($indvb1 eq $femelleb1) or ($indvb2 eq $femelleb1))){
						#print OUT "ab\t";
					}
					elsif($classe eq "NN"){
						#print OUT "--\t";
					}
					else{
						#print OUT "NO\t";
					}
				}
				#print OUT "\n";
				#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxaa>\t";
				for(my $i = 4;$i<@tab;$i++){
					$tab[$i] =~ m/(\w+)\|(\w+)/;
					$indvb1 = $1;
					$indvb2 = $2;
					$classe = "$indvb1$indvb2";
					
					if((length($1) == 1) and (length($2) == 1)){
						my $chaine = "$indvb1$indvb2";
						my $GN = "$1$2";
						my @GN = split ('', $GN);
						my @ch = sort({$a cmp $b} @GN );
						my $GNT=join('', @ch);
						$classe = "$GNT";
					}
					
					if(($classe eq $femelle) and ($classe ne "NN")){
						#print OUT "aa\t";
					}
					elsif(($indvb1 ne $indvb2) and (($indvb1 eq $femelleb1) or ($indvb2 eq $femelleb1))){
						#print OUT "ab\t";
					}
					elsif($classe eq "NN"){
						#print OUT "--\t";
					}
					else{
						#print OUT "NO\t";
					}
				}
					print SUP "err:50	$line\n";
					$_50++;
			}
			# a0xab ou ab x a0 
			elsif($nbr_classes == 3){
				if($nbr_hetero == 1){
					#55# obs: AA x AA ,reel: err allele nul: A0 x AC ou AC x A0 => doubler
					if ($countHETERO <= ($countHOMO/2)){
						#print SUP "S: AAxAA_A0xAC	$line\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos\t<a0xab>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if($classe eq $femelle){
								#print OUT "a-\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
								#print OUT "b0\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq $femelleb1) or ($indvb2 eq $femelleb1))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						#print OUT "\n";
						#print OUT "$contig\t$pos\t$contig"."_"."$pos"."_i"."\t<abxa0>\t";
						for(my $i = 4;$i<@tab;$i++){
							$tab[$i] =~ m/(\w+)\|(\w+)/;
							$indvb1 = $1;
							$indvb2 = $2;
							$classe = "$indvb1$indvb2";
							
							if((length($1) == 1) and (length($2) == 1)){
								my $chaine = "$indvb1$indvb2";
								my $GN = "$1$2";
								my @GN = split ('', $GN);
								my @ch = sort({$a cmp $b} @GN );
								my $GNT=join('', @ch);
								$classe = "$GNT";
							}
							
							if($classe eq $femelle){
								#print OUT "a-\t";
							}
							elsif(($indvb1 eq $indvb2) and ($classe ne $femelle) and ($classe ne "NN")){
								#print OUT "b0\t";
							}
							elsif(($indvb1 ne $indvb2) and (($indvb1 eq $femelleb1) or ($indvb2 eq $femelleb1))){
								#print OUT "ab\t";
							}
							elsif($classe eq "NN"){
								#print OUT "--\t";
							}
							else {
								#print OUT "NO\t";
							}
						}
						print SUP "err:51	$line\n";
						$_51++;
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

%CLASSES = ();
%HOMO = ();
%HETERO = ();
%count = ();

@clas1 = ();
@clas2 = ();
@clasHOMO = ();
@clasHETERO = ();
@heterob1 = ();
@heterob2 = ();
@homob1 = ();
@homob2 = ();
@baseHetero = ();
@baseHetero1 = ();
%hachage = ();
%HOMO_count = ();

$nbr_classes = 0;
$countHOMO = 0;
$countHETERO = 0;
$nbr_hetero = 0;
$nbr_homo = 0;

}

print INFO "indiv_marker = $indiv_marker\n";
print INFO "group_marker = $group_marker\n";
my $total = $indiv_marker+$group_marker;
print INFO "Nombre total de marqueurs traités : $total\n";

print INFO "1 = "." $_1 ;"; print INFO "2 = "." $_2 ;"; print INFO "3 = "." $_3 ;"; print INFO "4 = "." $_4 ;"; print INFO "5 = "." $_5 ;"; print INFO "6 = "." $_6 ;"; print INFO "7 = "." $_7 ;"; print INFO "8 = "." $_8 ;"; print INFO "9 = "." $_9 ;"; print INFO "10 = "." $_10 ;"; print INFO "11 = "." $_11 ;"; print INFO "12 = "." $_12 ;"; print INFO "13 = "." $_13 ;"; print INFO "14 = "." $_14 ;"; print INFO "15 = "." $_15 ;"; print INFO "16 = "." $_16 ;"; print INFO "17 = "." $_17 ;"; print INFO "18 = "." $_18 ;"; print INFO "19 = "." $_19 ;"; print INFO "20 = "." $_20 ;"; print INFO "21 = "." $_21 ;"; print INFO "22 = "." $_22 ;"; print INFO "23 = "." $_23 ;"; print INFO "24 = "." $_25 ;"; print INFO "25 = "." $_25 ;"; print INFO "26 = "." $_26 ;"; print INFO "27 = "." $_27 ;"; print INFO "28 = "." $_28 ;"; print INFO "29 = "." $_29 ;"; print INFO "30 = "." $_30 ;"; print INFO "31 = "." $_31 ;"; print INFO "32 = "." $_32 ;"; print INFO "33 = "." $_33 ;"; print INFO "34 = "." $_34 ;"; print INFO "35 = "." $_35 ;"; print INFO "36 = "." $_36 ;"; print INFO "37 = "." $_37 ;"; print INFO "38 = "." $_38 ;"; print INFO "39 = "." $_39 ;"; print INFO "40 = "." $_40 ;"; print INFO "41 = "." $_41 ;"; print INFO "42 = "." $_42 ;"; print INFO "43 = "." $_43 ;"; print INFO "44 = "." $_44 ;"; print INFO "45 = "." $_45 ;"; print INFO "46 = "." $_46 ;"; print INFO "47 = "." $_47 ;"; print INFO "48 = "." $_48 ;"; print INFO "49 = "." $_49 ;"; print INFO "50 = "." $_50;"; print INFO "51 = "." $_51 ;";print INFO "52 = "." $_52 ;"; print INFO "53 = "." $_53 ;"; print INFO "54 = "." $_54;"; print INFO "55 = "." $_55 ;";

close F;
close OUT;
close SUP;
close INFO;



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
  the Free Software Foundation; either version 3 of the License, or
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
