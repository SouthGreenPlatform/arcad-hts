#!/bin/env perl

#  
#  Copyright 2014 INRA-CIRAD-UM2
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

cleanMapping.pl - Remove reads that have not correctly mapped

=head1 DESCRIPTION

Work for bwa outfiles

Remove reads that have not correctly mapped, including :
1) soft-clipped
2) number of mismatches > maxMismatch, 
3) number of indels > maxIndels
4) all indels larger than maxIndelsSize

Finally, the user can choose if whether or not the reads that mapped at multiple places should be excluded :
if multipleMap = "true", the reads are kept and if multipleMap = "false", the reads are removed.

If reads are pair-end, if one pair is removed the other is also removed even if it doesn't fullfill the criteria above.

Create a {input sam}_reads_removed.txt containing the name of the removed reads.

There is two solutions for executing this script, users can use either a single bam file (--bam) in input or a folder 
and a pattern (*.bam) to specified a group of BAM files. Output files have same name than input BAM with extension .clean

=head1 SYNOPSIS / USAGE

cleanMapping.pl --bam algnSam [--maxMismatch maxMismatch] [--maxIndels maxIndels] [--maxIndelsSize maxIndelsSize] [--multipleMap multipleMap]

cleanMapping.pl --folder folder_algnSam/ [--maxMismatch maxMismatch] [--maxIndels maxIndels] [--maxIndelsSize maxIndelsSize] [--multipleMap multipleMap]

cleanMapping.pl --folder folder_algnSam/ --pattern *group1.sam [--maxMismatch maxMismatch] [--maxIndels maxIndels] [--maxIndelsSize maxIndelsSize] [--multipleMap multipleMap]

=cut

my $version="";

use strict;
use warnings;
use Carp qw (cluck confess croak);
use Pod::Usage;
use Getopt::Long;
use File::Path;
use File::Basename;
use File::Spec::Functions;

use lib '/NAS/arcad_data/Softs/tags';
use Modules::Files::Files;
use Modules::Config::Softwares;

sub help
{ 
	my $msg = shift;
	
	print "-----------------------------------------------------------------------------\n"
	."-- ERROR MESSAGE --\n"."$msg\n".
	"-----------------------------------------------------------------------------\n"
	if( defined  $msg );
	pod2usage(0);
}

help "Not enough arguments" unless (@ARGV >= 1);

=head1 OPTIONS

=head2 required parameters

=over 2

=item B<--bam> (bam file or folder) 

BAM folder in which to find all BAM files.

=back

=head2 optional arguments

=over

=item B<[--maxMismatch]> ([treshold of Mismatch]): DEFAULT 5

Remove read that number of mismatch is higher than this.

=item B<[--multipleMap]> ([Exclude reads that mapped at multiple places <true/false>]): DEFAULT false

this option is a boolean, if option --multiMap is specified in arguments, value will be true.
Otherwise value is set to false. If 'false' reads that mapped on several places are remove, and not if 'true'.

=item B<[--maxIndels]> ([maximum number of allowed Indel]): DEFAULT 2

=item B<[--maxIndelsSize]> ([Treshold for Indel size]): DEFAULT 4


=item B<[--pattern]> ([output folder]): DEFAULT *.sam

This option is used to liste all SAM files to use for the cleaning. 
This option is effective only if parameter --sam is a folder.

=item B<[--subdirectory|--nosubdirectory]> 

--subdirectory   when users want SAM files in specified directory and it's sub directories.
--nosubdirectory when users want only SAM files in specified directory.
This option is effective only if parameter --sam is a folder.

=back

=cut

my $bami;
my $maxMismatch = 5;
my $maxIndels = 2;
my $maxIndelsSize = 4;
my $multipleMap;
my $out = '.';

my $pattern = '*.bam';
my $subdirectory = '';
my $ra_files;

GetOptions(
	"bam=s"           => \$bami,
	"maxMismatch=i"   => \$maxMismatch,
	"maxIndels=i"     => \$maxIndels,
	"maxIndelsSize=i" => \$maxIndelsSize,
	"multipleMap"     => \$multipleMap,
	"output|o=s"      => \$out,
	
	"pattern=s"       => \$pattern,
	"subdirectory!"   => \$subdirectory,
	"help|?|h" => sub{ help(); }
);

if( !-e $bami ){
	help "$! $bami";
}
elsif( -f $bami )
{
	push @$ra_files, $bami;
}
elsif( -d $bami )
{
	$ra_files = Files->getFiles( $bami, $pattern, $subdirectory );
	help "No BAM file(s) in input" if( 0 >= scalar @$ra_files );
}
else{
	help "$! $bami";
}

#-----------------------------------------------------
# CREATION DOSSIER TMP_CLEAN/
#-----------------------------------------------------
mkdir $out if(defined($out) && !(-e $out));
my $outDir = $out;
$outDir = $outDir . "/tmp_clean/";
mkdir ($outDir) or die ("Erreur creation repertoire\n");
#-----------------------------------------------------
# EXECUTABLES path
#-----------------------------------------------------
my $SAMTOOLS_PATH = &$Softwares::SAMTOOLS_PATH;
#-----------------------------------------------------
# BAM TO SAM : BAM ds TMP_CLEAN/
#-----------------------------------------------------
foreach my $bam (@$ra_files)
{
	my $sam = $bam;
	$sam =~ s/\.bam$/\.sam/;	
	my $rh_qsub = &$Softwares::RH_QSUB;
	$$rh_qsub{N} = 'samtools_view';
	my $qsub = &Softwares::make_qsub_command($rh_qsub);
	$out = catfile($outDir, basename $sam);
	$? = system( "$qsub '$SAMTOOLS_PATH view -h $bam > $out'" );
	throw Error::Simple("Job 'samtools_view' failed, exit code $?") if( !-e $sam && $? != 0 );
#-----------------------------------------------------

	my $clean = catfile($outDir, basename $sam . '.clean');

	## if pair-end data rm the pair of a badly mapped reads
	my %RM_READS;

	open F, "$out" or confess("$!");
	my $cn = 0;
	while( my $line = <F>)
	{
		chomp($line);

		next if($line =~ /^\@SQ/ || $line =~ /^\@PG/ || $line =~ /^\@HD/ || $line =~ /^\@RG/);

		my @a = split(" ", $line);
		my $cigar = $a[5];
		my $NREF = $a[6];
		my $name = $a[0];
		
		if($cigar eq "*")
		{
			if(!exists($RM_READS{$name})){
				$RM_READS{$name} = 1;
			}
			next;
		}
		
		## Map at multiple place or two pairs don't mach in the same reference seq
		if(undef $multipleMap || !$multipleMap)
		{
			if($line !~ /X0\:i\:1/)
			{
				if(!exists($RM_READS{$name})){
					$RM_READS{$name} = 1;
				}
				next;
			}
			if($NREF ne "=" && $NREF ne "*")
			{
				if(!exists($RM_READS{$name})){
					$RM_READS{$name} = 1;
				}
				next;
			}
		}

		### check the number of mismatches
		my $mismatch = "";
		foreach my $it (@a)
		{
			if($it =~ /^NM\:/){
				$mismatch = $it;
			}
		}
		
		help ("$line\n") if($mismatch eq "");
		
		my @mis = split (":", $mismatch);
		if($mis[2] > $maxMismatch)
		{
			if(!exists($RM_READS{$name})){
				$RM_READS{$name} = 1;
			}
			next;
		}

		### check the cigar string
		## rm soft clipped reads or other strange thing in the cigar
		if($cigar =~ /H/ || $cigar =~ /N/ || $cigar =~ /S/)
		{
			if(!exists($RM_READS{$name})){
				$RM_READS{$name} = 1;
			}
			next;
		}

		my @cig = split(//, $cigar);
		my @ind;
		
		## retrieve effectif of each types
		foreach my $op (@cig)
		{
			if($op !~ /[0-9]/){ 
				push(@ind, $op);
			}
		}

		## split the cigar string
		@cig = split(/M|S|I|D/, $cigar);

		my $numberInDels = 0;
		my $sizeInDels = 0;
		my $c = 0;
		
		foreach my $op (@cig)
		{
			if($ind[$c] eq "I")
			{
				$numberInDels = $numberInDels + 1;
				for(my $i = 0; $i<$op; $i++){
					$sizeInDels = $sizeInDels + 1;
				}
			}
			if($ind[$c] eq "D")
			{
				$numberInDels = $numberInDels + 1;
				for(my $i = 0; $i<$op; $i++){
					$sizeInDels = $sizeInDels + 1;
				}
			}
			$c = $c + 1;
		}

		if($numberInDels > $maxIndels || $sizeInDels > $maxIndelsSize)
		{
			if(!exists($RM_READS{$name})){
				$RM_READS{$name} = 1;
			}
		}
	}
	close F;

	open OUT,">$clean" or confess "$!\n$clean\n";

	open F, "$out";
	while( my $line = <F>)
	{
		chomp($line);
		
		if($line =~ /^\@PG/ || $line =~ /^\@SQ/ || $line =~ /^\@HD/ || $line =~ /^\@RG/)
		{
			print OUT "$line\n";
			next;
		}

		my @a = split(" ", $line);
		my $name = $a[0];
		if(!exists($RM_READS{$name})){
			print OUT "$line\n";
		}
	}
	close F;
	close OUT;

	#-----------------------------------------------------
	# SAM TO BAM
	#-----------------------------------------------------
	my $outb = $bam;
	$outb =~ s/\.bam$/\.bam.clean/;	
	#my $rh_qsub = &$Softwares::RH_QSUB;
	#my $qsub = &Softwares::make_qsub_command($rh_qsub);
	#$$rh_qsub{N} = 'samtools_view';	
	$? = system( "$qsub '$SAMTOOLS_PATH view -bS $clean > $outb'" );
	throw Error::Simple("Job 'samtools_view' failed, exit code $?") if( !-e $outb && $? != 0 );
	#-----------------------------------------------------

	open O, '>', $sam . '_reads_removed.txt' or confess "$!\n";
	foreach my $n (keys %RM_READS){
		print O "$n\n";
	}
	close O;

	%RM_READS = ();
}

#-----------------------------------------------------
# EFFACEMENT DE TMP_CLEAN
#-----------------------------------------------------
if (-e $outDir) {
	rmtree(["$outDir"]);
}

exit 0;
	
=pod

=head1 AUTHORS

Created by Benoit Nabholz (benoit.nabholz@gmail.com)
Modified by Felix Homa & Gautier Sarah

=head1 VERSION

Version 1.0.1

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

=cut
