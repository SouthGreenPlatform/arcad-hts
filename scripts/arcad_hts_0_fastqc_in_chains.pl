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

=head1 FastQC in chains


Launch FastQC on all fastq (and fastq.gz if specified) files in a directory and its arborescence

=head1 SYNOPSIS

fastqc_in_chains.pl  -i input_folder [-o output_folder] [-q queue][-gz|-nogz] [-sub|-nosub] [-extract|-noextract]

=cut

use strict;
use Carp qw (confess);
use warnings;
use Readonly;
use Getopt::Long;
use Pod::Usage;

use lib '/NAS/arcad_data/Softs/tags';
use Modules::Config::Softwares;

=pod

=head1 OPTIONS


=head2 Parameters

=over 4

=item B<[-i]> ([input_folder]): 

The directory containing the fastq files in it or in its subdirectory


=item B<[-q]> ([queue]): 

The queue to launch FastQC. 

=head3 Default arcad.q

=item B<[-o]> ([output_folder]): 
If provided, FastQC outputs will be gather in this directory. 

=head3 Default behavior is to have the FastQC output in the same directory as the fastq file


=item B<[-gz|-nogz]> ([flag to process also fastq.gz file])
-gz will make the script to launch FastQC also on the fastq.gz files.
-nogz will make the script to launch FastQC only on the fastq files

=head3 Default: -gz is on


=item B<[-sub|nosub]> ([flag to check subdirectory])
-sub will allow to launch FastQC on fastq files present in subdirectories
-nosub disable this possibility

=head3 Default: -sub is on


=item B<[-extract|noextract]> ([flag to set if data must be keeped compressed or not])
-extract results will be extracted
-noextracted results will be kept as zipped files

=head3 Default:extract is on


=cut

my $version="";

unless (@ARGV)
{
	pod2usage({
		-message => 'No option found',
		-verbose => 0,
		});
	exit();
}

my $FASTQC_PATH = &$Softwares::FASTQC_PATH;

#options processing
my ($man, $help, $debug, $input, $output, $gz, $subdirectory,$extract, $queue);
$gz=1;
$subdirectory=1;
$extract=1;
$queue="normal.q";

# parse options and print usage if there is a syntax error.
#--- see doc at http://perldoc.perl.org/Getopt/Long.html
GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "input|i=s"  => \$input,
           "output|o=s" => \$output,
           "gz!"        => \$gz,
           "sub!"       => \$subdirectory,
           "extract!"   => \$extract,
           "queue|q=s"  => \$queue
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

#if sub is on we looked for all fastq files in the arborescence
my $list_command;
if($subdirectory)
{
  $list_command="find -L $input -type f -name \"*.fastq\" ";
  #if gz is on, we also looked for fastq.gz file
  if($gz)
  {
    $list_command=$list_command."; find -L $input -type f -name \"*.fastq.gz\"";
  }
}
else
{
  $list_command="ls -1 $input/*.fastq ";
  if($gz)
  {
    $list_command=$list_command."; ls -1 $input/*.fastq.gz";
  }
}

my @file_list = `$list_command`;
chomp(@file_list);
my $nbfile=scalar @file_list;
my $thread=1;

if(defined($debug))
{
  print "Liste of files : ", join("\n", @file_list), "\n";
  print "$nbfile files are ready to be processed\n";
  print "COMMAND USED : $list_command\n";
}

if($nbfile>11)
{
  $thread=12;
}
else
{
  $thread=$nbfile;
  chomp $thread;
}

my $option="";
if (defined($output))
{
  $option=$option." -o $output";
}
if(!$extract)
{
  $option=$option." --noextract";
}

if(defined($debug))
{
  print "FastQC command\n";
  print("qsub -q $queue -b yes -V -cwd -N FastQC_in_chains -pe parallel_fill $thread '$FASTQC_PATH -t $thread $option @file_list'\n");
}

system("qsub -q $queue -b yes -V -cwd -N FastQC_in_chains -pe parallel_fill $thread '$FASTQC_PATH -t $thread $option @file_list'");

=pod

=head1 AUTHORS

Gautier SARAH (INRA), gautier.sarah@supagro.inra.fr

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
 
