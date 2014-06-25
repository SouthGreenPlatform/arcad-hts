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

eliminer_distortion_kideux.pl - eliminates distorted SNPs with probability below 0.01 
and a probability less than 0.05 (for a window of 10 markers)

=head1 DESCRIPTION

Eliminates distorted SNPs with probability below 0.01 
and a probability less than 0.05 (for a window of 10 markers)

=head1 SYNOPSIS / USAGE

eliminer_distortion_kideux.pl	[--help] [--man] [--version]
[-i input-file] \
[-o output] \


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
input file, (output of kideux.pl script)

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

my (@tab, $marker, $dest, $segregation, @chr, @pos, @line, @head);
my $e = ""; my $ee = ""; my $eee = "";
my $head = <F>;

print OUT "$head";
my $d = 0;
my $f = 0;
my $i = 0;

while(my $line =<F>){
	chomp $line;
	@tab = split (/\t/, $line);
	@head = split (/\t/, $head);
	push (@chr, $tab[0]);

	$marker = $tab[2];
	$segregation = $tab[3];
	
	$e = $tab[-3];
	$ee = $tab[-2];
	$eee = $tab[-1];
	
	if($chr[$i] eq $chr[$i-1]){
		$f++;
		if($f< 10){
			if($eee eq "*"){
				#print "$line\n";
				$d++;
			}
			elsif(($ee eq "*") and ($eee ne "*")){
				#print "$line\n";
				$d++;
			}
			elsif(($e eq "*") and ($ee ne "*") and ($eee ne "*")){
				$d++;
				push(@line, $line);
			}
			else{
				print OUT "$line\n";
			}
			
		}elsif( $f == 10){
			if($eee eq "*"){
				#print "$line\n";
				$d++;
			}
			elsif(($ee eq "*") and ($eee ne "*")){
				#print "$line\n";
				$d++;
			}
			elsif(($e eq "*") and ($ee ne "*") and ($eee ne "*")){
				$d++;
				push(@line, $line);
			}
			else{
				print OUT "$line\n";
			}
			if($d >= 2){
				foreach my $val(@line){
					print OUT "$val\n";
				}
			}
			else{
				foreach my $val(@line){
					print "$val\n";
				}
			}
			#print OUT "10\n";
			$d = 0;
			$f = 0;
			@line = ();
		}
		else{
			#print OUT "PLUS DE 10\n";
		}
	}
	#Un autre chromosome
	else{
		if($d >= 2){
			foreach my $val(@line){
				print OUT "$val\n";
			}
		}else{
			foreach my $val(@line){
				print "$val\n";
			}
		}
		$f = 0;
		$d = 0;
		@line = ();
		$f++;
		#print OUT "Chromosome\n";
		if($eee eq "*"){
			#print "$line\n";
			$d++;
		}
		elsif(($ee eq "*") and ($eee ne "*")){
			#print "$line\n";
			$d++;
		}
		elsif(($e eq "*") and ($ee ne "*") and ($eee ne "*")){
			push(@line, $line);
			$d++;
		}
		else{
			print OUT "$line\n";
		}
	}
	$i++;
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
