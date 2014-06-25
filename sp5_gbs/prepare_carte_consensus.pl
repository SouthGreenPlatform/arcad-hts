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

prepare_carte_consensus.pl - prepare input for consensus map. Coding compatible with Carthagene tool

=head1 DESCRIPTION

Prepare input for consensus map.Adds a column "phase" in the input file (abcd coding-file)

Coding compatible with Carthagene tool


=head1 SYNOPSIS / USAGE

prepare_carte_consensus.pl	[--help] [--man] [--version] 

[-i tabulated-file (Carthagene output for grouping)] \

[-f abcd coding-file] \

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
input file (Carthagene output for grouping)


=item B<[-f]> ([input_file]): 
abcd coding-file


=item B<[-o]> ([prefix_output]): 
prefix_output 


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

my ($man, $help, $version, $input1, $input2, $prefix_output);
$man = 0;
$help = 0;
$version = 0;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"input|i=s"   => \$input1,
	"file|f=s"   => \$input2,
	"output|o=s"  => \$prefix_output	
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
if ($version) {print "1.0.0\n";exit}

my $output = "$prefix_output"."_consensus.txt";

if(defined($input1) && ! -e $input1 )
{
	warn "$input1 not exists ";
	exit(0);
	
}
if(defined($input2) && ! -e $input2 )
{
	warn "$input2 not exists ";
	exit(0);
	
}
if(defined($output) && -e $output )
{
	warn "$output already exists ";
	exit(0);
	
}

if(defined($output) && -e $output )
{
	warn "$output already exists ";
	exit(0);
	
}

open (F1, $input1)or die "cannot open the file $input1\n";
open (F2, $input2)or die "cannot open the file $input2\n";
open(OUT, ">$output") or die"cannot open the file $output\n";

my (@markerM, @markerF, $femele, $male, $head);

while(my $line =<F1>){
	chomp $line;
	my @tab1 = split (/\t/, $line);
	my $parent = $tab1[0];
	my $marker = $tab1[1];
	
	if ($parent eq "male"){
		push (@markerM, $marker);
	}
	elsif($parent eq "female"){
		push (@markerF, $marker);
	}
}

close F1;
$head = <F2>;
my @head = split (/\t/, $head);
	for(my $i = 0;$i<3;$i++){
		print OUT "$head[$i]\t";
	}
	print OUT "Phasing\t";
	for(my $i = 3;$i<@head;$i++){
		print OUT "$head[$i]\t";
	}
while(my $line =<F2>){
	chomp $line;
	my @tab1 = split (/\t/, $line);
	my $marker1 = $tab1[2];
	my $miror = $marker1."_m";
	#foreach my $val (@markerF){
		
			#if($val eq "$marker1"){
				#$femele = "0";
			
			#}
			#elsif($val eq "$marker1"."_m"){
				
				#$femele = "1";
			#}
			#elsif(($val ne "$marker1"."_m") and ($val ne "$marker1")){
				#next;
			#}
			
			if (grep $_ eq $marker1, @markerF) {
				$femele = "0"; 
			}
			elsif(grep $_ eq $miror, @markerF){
				$femele = "1"; 
			}
			else{
				$femele = "-"; 
			}
			#$femele = "0" if grep {$_ eq $marker1} @markerF;
			#$femele = "1" if grep {$_ eq $miror} @markerF;
			#$femele = "-" if grep {$_ ne $marker1} @markerF;
			#$femele = "-" if grep {$_ ne $miror} @markerF;
			#if (! exists $marker1 (@markerF)){
				
				
			#}
		#}
	#}

	#foreach my $val (@markerM){
		#if ( exists $marker1 (@markerM)){
			#if($val eq "$marker1"){
				
				#$male = "0";
			
			#}
			#elsif($val eq "$marker1"."_m"){
				
				#$male = "1";
			#}
			#elsif(($val ne "$marker1"."_m") and ($val ne "$marker1")){
				#next;
			#}
			#else{
				##$male = "-";
			#}
			if (grep $_ eq $marker1, @markerM) {
				$male = "0"; 
			}
			elsif(grep $_ eq $miror, @markerM){
				$male = "1"; 
			}
			else{
				$male = "-"; 
			}
			#$male = "0" if grep {$_ eq $marker1} @markerM;
			#$male = "1" if grep {$_ eq $miror} @markerM;
			#$male = "-" if grep {$_ ne $marker1} @markerM;
			#$male = "-" if grep {$_ ne $miror} @markerM;
		#}
	#}
	if(($femele ne '') and ($male ne '')){
		for(my $i = 0;$i<3;$i++){
			print OUT "$tab1[$i]\t";
		}
		print OUT "{"."$femele"."$male}\t";
		for(my $i = 3;$i<@tab1;$i++){
			print OUT "$tab1[$i]\t";
		}
		print OUT "\n";
	}
	else{
		next;
	}
	
	$femele = '';
	$male = '';

}

close F2;
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
