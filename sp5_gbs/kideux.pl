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

Kideux.pl - do a chi-square test for markers 

=head1 DESCRIPTION

Do a chi-square test for markers 

=head1 SYNOPSIS / USAGE

Kideux.pl	[--help] [--man] [--version]
[-i input-file] \
[-o output] \


example of input file :

Chromosome	POS	Marker	Segregation	Child1	Child2	Child3

chr1	192497	chr1_192497	<abxaa>	aa	aa	aa

chr1	235611	chr1_235611	<abxaa>	aa	--	aa

chr1	235683	chr1_235683	<abxaa>	aa	--	aa


=cut

use strict;
use Pod::Usage;
use Getopt::Long;
use Carp qw (cluck confess croak);
use Statistics::Distributions ; 

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

=over 5

=item B<[-i]> ([input_file]): 
input file


=item B<[-o]> ([output]): 
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

if(defined($input) && ! -e $input )
{
	warn "$input not exists ";
	exit(0);
	
}
if(not defined($output))
{
	warn "you forgot the -o parameter ";
	exit(0);
	
}
if(defined($output) && -e $output )
{
	warn "$output already exists ";
	exit(0);
	
}

open (F, $input) or die "cannot open the file $input\n";
open(OUT, ">$output") or die"cannot open the file $output\n";

my (@tab,$NO, $e1, $e2, $e3, $contig, $pos, $base_ref, $base_alt, $femelle, $male, $classe, $segregation, $total, $aa, $ab, $bb, $amoins, $bmoins, $ac, $ad, $bc, $bd, $NN, $bzero, $chi2, $ddl, $proba, $DP);
$total = $aa = $ab = $bb = $amoins = $bmoins = $ac = $ad = $bc = $bd = $NN = $bzero = $chi2 = $ddl = $proba = $DP = $NO = 0;
$e1 = ""; $e2 = ""; $e3 = "";

my $head = <F>;
$head =~ s/\n//g;
print OUT $head;
print OUT "	aa	ab	bb	a-	b-	ac	ad	bc	bd	--	b0	total	DP	chi2	ddl	proba	e	ee	eee\n";
while(my $line =<F>){
	chomp $line;
	@tab = split (/\t/, $line);
	$contig = $tab[0];
	$pos = $tab[1];
	#$segregation = $tab[5];
	$segregation = $tab[3];
	#Output format :
	#Indiv		aa	ab	bb	a-	b-	ac	ad	bc	bd	--	b0	total	DA	chi2	ddl	proba	e	ee	eee
	
				#counting each class
				for(my $i = 3;$i<@tab;$i++){
					$classe = $tab[$i];
					if($classe eq "aa"){
						$aa++;
					}
					elsif($classe eq "ab"){
						$ab++;
					}
					elsif($classe eq "bb"){
						$bb++;
					}
					elsif($classe eq "a-"){
						$amoins++;
					}
					elsif($classe eq "b-"){
						$bmoins++;
					}
					elsif($classe eq "ac"){
						$ac++;
					}
					elsif($classe eq "ad"){
						$ad++;
					}
					elsif($classe eq "bc"){
						$bc++;
					}
					elsif($classe eq "bd"){
						$bd++;
					}
					elsif($classe eq "--"){
						$NN++;
					}
					elsif($classe eq "b0"){
						$bzero++;
					}
					else{
						$NO++;
					}
				}
				$total = $aa + $ab + $bb + $amoins + $bmoins + $ac + $ad + $bc + $bd + $NN + $bzero ; 
				$DP = $total - $NN ;
				if(($segregation eq "<abxaa>") or ($segregation eq "<aaxab>")){
					$ddl = 1;
					$chi2 = 2*(($aa-$DP*0.5)**2)/$DP/0.5;
					$proba = Statistics::Distributions::chisqrprob($ddl, $chi2);
				}elsif(($segregation eq "<a0xab>") or ($segregation eq "<abxa0>")){
					$ddl = 2;
					$chi2 = ($ab-$DP*0.25)**2/$DP/0.25+($bzero-$DP*0.25)**2/$DP/0.25+($amoins-$DP*0.5)**2/$DP/0.5;
					$proba = Statistics::Distributions::chisqrprob($ddl, $chi2);
				}
				elsif($segregation eq "<abxab>"){
					$ddl = 2;
					$chi2 = ($aa-$DP*0.25)**2/$DP/0.25+($bb-$DP*0.25)**2/$DP/0.25+($ab-$DP*0.5)**2/$DP/0.5;
					$proba = Statistics::Distributions::chisqrprob($ddl, $chi2);
				}
				elsif($segregation eq "<abxac>"){
					$ddl = 3;
					$chi2 = ($aa-$DP*0.25)**2/$DP/0.25+($ac-$DP*0.25)**2/$DP/0.25+($ab-$DP*0.25)**2/$DP/0.25+($bc-$DP*0.25)**2/$DP/0.25;
					$proba = Statistics::Distributions::chisqrprob($ddl, $chi2);
				}
				elsif($segregation eq "<abxcd>"){
					$ddl = 3;
					$chi2 = ($ac-$DP*0.25)**2/$DP/0.25+($ad-$DP*0.25)**2/$DP/0.25+($bc-$DP*0.25)**2/$DP/0.25+($bd-$DP*0.25)**2/$DP/0.25;
					$proba = Statistics::Distributions::chisqrprob($ddl, $chi2);
				}
				else{
					#just for test
					$ddl = "NO";
					$chi2 = "NO";
					#$proba = pchisq();
				}
				#put a star (*) for each interval probability
				if($proba < 0.05){
					$e1 = "*";
				}
				else{
					$e1 = 0;
				}
				if($proba < 0.01){
					$e2 = "*"
				}
				else{
					$e2 = 0;
				}
				if($proba < 0.001){
					$e3 = "*";
				}
				else{
					$e3 = 0;
				}
				print OUT "$line	$aa	$ab	$bb	$amoins	$bmoins	$ac	$ad	$bc	$bd	$NN	$bzero	$total	$DP	$chi2	$ddl	$proba	$e1	$e2	$e3\n";

$total = $aa = $ab = $bb = $amoins = $bmoins = $ac = $ad = $bc = $bd = $NN = $bzero = $chi2 = $ddl = $proba = $DP = 0;
$e1 = ""; $e2 = ""; $e3 = "";

}

close F;
close OUT;


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
