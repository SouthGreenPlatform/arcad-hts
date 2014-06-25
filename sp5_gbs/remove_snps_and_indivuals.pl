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

remove_snps_and_indivuals.pl - removes individuals and SNPs with more than X missing data 

=head1 DESCRIPTION

Removes individuals and SNPs with more than X missing data 

Generates a file "prefix_removed_indiv.csv" with a list of eliminated individuals

=head1 SYNOPSIS / USAGE

remove_snps_and_indivuals.pl	[--help] [--man] [--version]

[-i csv tabulated-file] \
[-o prefix_output] \
[-a maximum missing data per SNP] \
[-b maximum missing data per individual] \


=cut

use strict;
use warnings;
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
csv file

=item B<[-o]> ([prefix_output]): 
prefix for output file

=item B<[-a]> ([missing data per SNP]): 
maximum missing data per SNP
=cut

=item B<[-b]> ([missing data per individual]): 
maximum missing data per individual


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

my ($man, $help, $version, $input, $prefix_output, $NA_SNP, $NA_indiv);
$man = 0;
$help = 0;
$version = 0;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"input|i=s"   => \$input,
	"a|a=s"   => \$NA_SNP,
	"b|b=s"   => \$NA_indiv,
	"output|o=s"  => \$prefix_output	
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
if ($version) {print "1.0.0\n"; exit}

if($NA_SNP eq '')
{
	print "You have not specified parameter -a\n";
	exit;
}
if($NA_indiv eq '')
{
	print "You have not specified parameter -b\n";
	exit;
}

my $out = $prefix_output."_tmp.csv";
my $sup = $prefix_output."_removed_indiv.csv";
if(! -e $input){ print "the file $input does not exist\n";}

open (F, $input) or die "cannot open the file $input\n";
open(OUT, ">$out") or die"cannot open the file $out\n";
open(SUP, ">$sup") or die"cannot open the file $sup\n";


my $head=<F>;
chomp $head;
print OUT "$head\n";
my $line;
while($line = <F>){
	chomp $line;
	my @tab = split (/\t/, $line);
	next if eof;
	if ($tab[$#tab] <= $NA_SNP){
		for(my $i = 0;$i<$#tab;$i++){
			print OUT "$tab[$i]\t";
		}
		print OUT "\n";
		
	}
}

close F;
close OUT;

my (@tab,$contig, $pos, $base_ref, $base_alt, $NA, %hash);
open(OUT, "+<$out") or die"cannot open the file $out\n";

while(my $line =<OUT>){
	chomp $line;
	@tab = split (/\t/, $line);
	for(my $i = 6;$i<@tab;$i++){
		if($tab[$i] eq "NN"){
			
			if ( defined $hash{$i} ){
				$hash{$i}++;
			}
			else{
				$hash{$i} = 1;
			}
		}
	}
}

print OUT "				0	0	";
for( my $j=6;$j<@tab;$j++){
    if ( defined $hash{$j} ){
        print OUT "$hash{$j}\t";
    }
    else{
        print OUT "0\t";
    }
}
close OUT;

my $out2 = $prefix_output."_tmp2.csv";
my $out3 = $prefix_output."_tmp3.csv";

open(OUT3, ">$out3") or die"cannot open the file $out3\n";

TransposerFichier("$out","$out2");

open(OUT2, "$out2") or die"cannot open the file $out2\n";
my $i = 0;
while(my $line =<OUT2>){
	chomp $line;
		@tab = split (/\t/, $line);
		if($i >=4 ){
			if ($tab[$#tab] <= $NA_indiv){ 
				for(my $i = 0;$i<$#tab;$i++){
					print OUT3 "$tab[$i]\t";
				}
				print OUT3 "\n";
			}
			else{
				print SUP "$tab[0]=> $tab[$#tab]\n";
			}
		}
		else{
			print OUT3 "$line\n";
		}
		$i++;
	}

my $out4 = $prefix_output.".csv";

close OUT;
close OUT2;
close OUT3;

TransposerFichier("$out3","$out4");

system ("rm $out");
system ("rm $out2");
system ("rm $out3");


#============================================
# TransposerFichier
# Transpose tabulated file
#convert columns in lines
#============================================
sub TransposerFichier {
  my ( $FichierTabuleOriginal, $FichierTranspose ) = @_;
 
  my %HashTranspose;
 
  # Reading tabulated file
  open( my $FH, '<', $FichierTabuleOriginal ) 
    or die("Can not open the file $FichierTabuleOriginal\n");
 
  while ( my $Line = <$FH> ) {
    chomp $Line;
    my @data = split( /\t/, $Line );
    for ( my $i = 0; $i < scalar(@data); $i++ ) {
      $HashTranspose{$i} .= $data[$i] . "\t";
    }
  }
  close($FH);
 
  # Creating transposed file
  open( my $FHTranpose, '>', $FichierTranspose ) 
    or die("Can not create the file $FichierTranspose\n");
 
  foreach ( sort { $a <=> $b } keys %HashTranspose ) {
    $HashTranspose{$_} =~ s{\t$}{};
    print {$FHTranpose} "$HashTranspose{$_}\n";
  }
  close($FHTranpose);
 
  return $FichierTranspose;
}

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
