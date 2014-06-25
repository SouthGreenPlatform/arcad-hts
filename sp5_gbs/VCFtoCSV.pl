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

VCFtoCSV.pl - transforms coding 01 of the genotypes (vcf file) into nucleotide coding ACGT (csv file)

=head1 DESCRIPTION

Transforms coding 01 of the genotypes into nucleotide coding and count the missing data per SNP (without parents) and per individual.

If a genotypic class is present only once, puts it into missing data

Generate a file "prefix_removed.csv"  containing removed SNPs 

=head1 SYNOPSIS / USAGE

VCFtoCSV.pl	[--help] [--man] [--version]
[-i VCF-file] \
[-o prefix_output] \


=cut

use strict;
use Pod::Usage;
use Getopt::Long;
use Carp qw (cluck confess croak);


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

VCF file


=item B<[-o]> ([prefix_output]): 

prefix for output file


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

my ($man, $help, $version, $input, $prefix_output);
$man = 0;
$help = 0;
$version = 0;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"input|i=s"   => \$input,
	"output|o=s"  => \$prefix_output	
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
if ($version) {print "1.0.0\n";exit}

my $out = $prefix_output.".csv";
my $sup = $prefix_output."_removed.csv";

if(! -e $input){ print "the file $input does not exist\n";}

open (F, $input) or die "cannot open the file $input\n";
open(OUT, ">$out") or die"cannot open the file $out\n";
open(SUP, ">$sup") or die"cannot open the file $sup\n";

my (@tab, @tab_base_alt, %CLASSES, %count, $classe, $nbr_classes, $chr, $pos, $base_ref, $base_alt, $alt_base1, $alt_base2, $alt_base3, $b1, $b2, $NN, %hash);

$NN =0;
$nbr_classes =0;
%count = ( );

print OUT "Chromosome	POS	REF	ALT	";
while(my $line = <F>){
	chomp $line;
	if($line =~ m/^#/){
		if ($line =~ m/^#CHROM/) {
		@tab = split (/\t/, $line);
			for(my $i = 9;$i<@tab;$i++){
				print OUT "$tab[$i]\t";
			}
			print OUT "\n";
		}
	}
	else{
		@tab = split (/\t/, $line);
		$chr = $tab[0];
		$pos = $tab[1];
		$base_ref = $tab[3];
		$base_alt = $tab[4];
		@tab_base_alt = split (',', $base_alt);
		$alt_base1 = $tab_base_alt[0];
		$alt_base2 = $tab_base_alt[1];
		$alt_base3 = $tab_base_alt[2];

		# keep only biallelic SNPs
		if(@tab_base_alt == 1){
			# count genotypic classes, i=9 => femele and i=10 => male; we start by i=11 => first child
			for(my $i = 11;$i<@tab;$i++){
				if ($tab[$i] =~ m/^(\d)\/(\d)*/){
					$classe = "$1$2";
					if(!exists($CLASSES{$classe})){
						$CLASSES{$classe} = 1;
						$nbr_classes++;
					}
					if(exists($CLASSES{$classe})){
						$count{$classe}++;
					}
				}
			}
			# keep SNPs that have more than one genotypic class in the offspring (remove the case: one genotypic class)
			if ($nbr_classes > 1){
				print OUT "$chr\t$pos\t$base_ref\t$base_alt\t";
				# write parent's genotypes
				for(my $i = 9;$i<11;$i++){
					if ($tab[$i] =~ m/^(\d)\/(\d)*/){
						if($1 == 0){
						$b1 = $base_ref;
						}
						if($1 == 1){
							$b1 = $alt_base1;
						}
						if($2 == 0){
							$b2 = $base_ref;
						}
						if($2 == 1){
							$b2 = $alt_base1;
						}
						my $chaine = "$b1$b2";
							# write genotypes in alphabetical order
							my @chaine = split ('', $chaine);
							my @ch = sort({$a cmp $b} @chaine );
							my $chano=join('', @ch);
							print OUT "$chano\t";
					}
					elsif($tab[$i] =~ m/^(\.)\/(\.)/){
						print OUT "NN\t";
						#$NN++;
						if ( defined $hash{$i} ){
							$hash{$i}++;
						}
						else{
							$hash{$i} = 1;
						}
					}
				}
				# write descendants's genotypes
				for(my $i = 11;$i<@tab;$i++){
					if ($tab[$i] =~ m/^(\d)\/(\d)*/){
						if($1 == 0){
						$b1 = $base_ref;
						}
						if($1 == 1){
							$b1 = $alt_base1;
						}
						if($2 == 0){
							$b2 = $base_ref;
						}
						if($2 == 1){
							$b2 = $alt_base1;
						}
						my $chaine = "$b1$b2";
						# If a genotypic class is present only once, transform it into NN (missing data)
						if($count{"$1$2"} == 1){
							$chaine = "NN";
							$NN++;
							$hash{$i}++;
						}
							# write genotypes in alphabetical order
							my @chaine = split ('', $chaine);
							my @ch = sort({$a cmp $b} @chaine );
							my $chano=join('', @ch);
							print OUT "$chano\t";
					}
					elsif ($tab[$i] =~ m/^(\.)\/(\.)/){
						print OUT "NN\t";
						$NN++;
						if ( defined $hash{$i} ){
							$hash{$i}++;
						}
						else{
							$hash{$i} = 1;
						}
					}
				}
				print OUT "$NN\n";
				$NN = 0;
			}
			else{
				print SUP "UneClasse	$line\n";
			}
		}
		else{
			print SUP "MultiAllélique	$line\n";
		}
	}

%CLASSES = ();
$nbr_classes = 0;
%count = ( );

}

print OUT "				";
for( my $j=9;$j<@tab;$j++){
    if ( defined $hash{$j} ){
        print OUT "$hash{$j}\t";
    }
    else{
        print OUT "0\t";
    }
}

close F;
close OUT;
close SUP;

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
