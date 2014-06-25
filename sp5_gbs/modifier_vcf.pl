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

modifier_vcf.pl - Modifies the VCF file to be accepted by PhaseByTransmission GATK module 

=head1 DESCRIPTION

Modifies the VCF file to be accepted by PhaseByTransmission GATK module 

Generates a pedigree file, will be used for the phasing

Generates à tabulated file, with list of added columns, will be used to remove this columns after

=head1 SYNOPSIS / USAGE

modifier_vcf.pl	[--help] [--man] [--version]
[-i VCF-file] \
[-o output] \


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


=item B<[-o]> ([output_file]): 

Output file (VCF)


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

my ($man, $help, $version, $input, $output, $ped, $indiv);
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
if ($version) {print "1.0.0\n";exit}


if(! -e $input){ print "the file $input does not exist\n";}

open (F, $input) or die "cannot open the file $input\n";
open(OUT, ">$output") or die"cannot open the file $output\n";
open(PED, ">$ped") or die "cannot open the file $ped\n";
open(INDIV, ">$indiv") or die "cannot open the file $indiv\n";


$ped = $output.".ped";
$indiv = $output."_Indiv_a_garder.txt";


while(my $line =<F>){
	chomp $line;
		if ($line =~ m/^##/) {
			print OUT "$line\n";
		}
		elsif($line =~ m/^#CHROM/ ){
			my @tab = split(/\t/, $line);
			for(my $i = 0;$i<12;$i++){
				print OUT "$tab[$i]\t";
				
			}
			print PED "FAM00Vd	$tab[9]	0	0	0	0\n";
			print PED "FAM00Vd $tab[10]	0	0	0	0\n";
			print PED "FAM00Vd	$tab[11]	$tab[10]	$tab[9]	0	0\n";
			print INDIV "$tab[9]\n$tab[10]\n$tab[11]\n";
			for(my $i = 12;$i<@tab;$i++){
				print OUT "Mo$i	Fa$i\t$tab[$i]	";
				print PED "FAM"."$i"."Vd	Mo$i	0	0	0	0\n";
				print PED "FAM"."$i"."Vd	Fa$i	0	0	0	0\n";
				print PED "FAM"."$i"."Vd	$tab[$i]	Fa$i	Mo$i	0	0\n";
				print INDIV "$tab[$i]\n";
			}
			print OUT "\n";
		}
		else{
			my @tab = split(/\t/, $line);
			for(my $i = 0;$i<12;$i++){
				print OUT "$tab[$i]\t";
				
			}
			for(my $i = 12;$i<@tab;$i++){
				print OUT "$tab[9]	$tab[10]	$tab[$i]\t";
			}
			print OUT "\n";
		}
	}

close F;
close OUT;
close PED;
close INDIV;


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
