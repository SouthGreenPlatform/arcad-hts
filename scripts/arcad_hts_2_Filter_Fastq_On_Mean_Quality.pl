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

Filter_Fastq_On_Mean_Quality - Filter FASTQ reads on mean quality

=head1 DESCRIPTION

Filter FASTQ reads on mean quality

=head1 SYNOPSIS / USAGE

Filter_Fastq_On_Mean_Quality.pl -f fastq_file -o output_file [-a ascii(default 33)] [-m mean_quality_threshold(default 30)] [-l length(default 35)]  

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Readonly;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Spec::Functions;

use lib '/NAS/arcad_data/Softs/tags';

use Modules::Files::Files;

sub help
{ 
	my $msg = shift;
	
	print "-----------------------------------------------------------------------------\n"
	."-- ERROR MESSAGE --\n"."$msg\n".
	"-----------------------------------------------------------------------------\n"
	if( defined  $msg );
	pod2usage(0);
}

pod2usage(0) unless (@ARGV);

=head1 OPTIONS

=head2 required parameters

=over 1

=item B<-f | --fastq> (Fastq file or directory)

Fastq file to proccess

=item B<-o | --output> (Output directory)

Directory in wich all output files will be sent.
Output files are same name as input with extension ".filtered"

=back

=head2 optional arguments

=over 1

=item B<-a | --ascii> DEFAULT 33

ascii code to determine the biginning of quality range

33 for Sanger

64 for Illumina 1.3 & 1.5

33 for Illumina 1.8 and higher

=item B<-m | --mean> DEFAULT 30

Minimum mean quality allowed for a read

=item B<-l | --length> DEFAULT 35

Minimum length allowed for a read

=item B<-? | --man | --help>

Print this help

=back

=cut

my $version="";

my ($man, $help, $ascii, $mean, $length, $fastq, $output) = (0, 0, 33, 30, 35, undef, './');
my $pattern = '*.fastq';
my $subdirectory = 0;

GetOptions("help|?"     => \$help,
           "man"        => \$man,
#           "debug:i"   => \$debug,
           "ascii|a=i"  => \$ascii, 
           "mean|m=i"   => \$mean,
           "length|l=i" => \$length,
           "fastq|f=s"  => \$fastq,
           "output|o=s" => \$output,
	   
	   "subdirectory!" => \$subdirectory,
	   "pattern|p=s" => \$pattern
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

help "No FASTQ file(s) found" unless( $fastq );

my $ra_files;

if( !-e $fastq ){
	help "$! $fastq";
}
elsif( -f $fastq )
{
	push @$ra_files, $fastq;
}
elsif( -d $fastq )
{
	$ra_files = Files->getFiles( $fastq, $pattern, $subdirectory );
	help "No FASTQ file(s) found" if( 0 >= scalar @$ra_files );
}
else{
	help "$! $fastq";
}

mkdir $output if(defined($output) && !(-e $output));
my $outputDir = $output;

foreach $fastq ( @$ra_files )
{
	$output = catfile($outputDir, basename $fastq . '.filtered');
	
	if ( open(my $file_handle, $fastq) && open(my $output_handle,">$output") )
	{
		my $score;
		while(<$file_handle>)
		{
			my $read=$_ . <$file_handle> . <$file_handle>; 
			$_ = <$file_handle>;
			chomp;
			map{ $score += ord($_)-$ascii} split(""); 
			print {$output_handle} "$read$_\n" if ($score/length($_)>= $mean && length($_)>=$length) ;
			$score= 0;
		}
		close($file_handle);
		close($output_handle);
	}
	else{
	     help "ERROR: failed to open file:\n$!\n";
	}
}

exit(0);

=head1 DEPENDENCIES

use Carp
use warnings
use Readonly
use Getopt::Long
use Pod::Usage

=head1 INCOMPATIBILITIES

<Fully compatble with all version of perl :)>

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 1

=item Gautier SARAH (INRA) gautier.sarah-at-cirad.fr

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

=cute)

=cut
