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

filtre_NO.pl - Eliminates the genotype markers taged"NO" by the "marker_to_ab.pl" script

=head1 DESCRIPTION

Eliminates the genotype markers taged "NO" by the "marker_to_ab.pl" script

these markers are due to recombination intra-markers or to phasing errors

=head1 SYNOPSIS / USAGE

filtre_NO.pl	[--help] [--man] [--version]
[-i input-file] \
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

my $sup = $output."_removed.txt";

open (F, $input)or die "cannot open the file $input\n";
open(OUT, ">$output") or die"cannot open the file $output\n";
open(SUP, ">$sup") or die"cannot open the file $sup\n";

my (@tab);
my $NO = 0;
while(my $line =<F>){
	chomp $line;
	@tab = split (/\t/, $line);
	#header
	if($line =~ m/^Chromosome/){
		@tab = split (/\t/, $line);
		print OUT "$line\n";
	}
	else{
		for(my $i = 4;$i<@tab;$i++){
			if($tab[$i] eq "NO"){
				$NO++;
			}
		}
		if($NO == 0){
			print OUT "$line\n";
		}
		else{
			print SUP "$line\n";
		}
	}
	$NO = 0;
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

=head1 TODO

Print warning message if input directory is same than output directory.
Because input data are overwritten.

=cut
