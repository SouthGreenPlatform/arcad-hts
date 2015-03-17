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

bamRealignerRecalibrate.pl - make (global) bam local realignment and base quality recalibration.

=head1 DESCRIPTION

According to GATK (2.0 verson) best practices, this script made LOCAL REALIGNMENT around indels in bam files.
Many BAM files can be treat at the same time, results can be sent in a unique bam or in several bams according 
to the number of in input.
In this script, we use GATK following modules :
RealignerTargetCreator
IndelRealigner
BaseRecalibrator
PrintReads

=head1 SYNOPSIS / USAGE

[% Not yet commented %]

=cut

use strict;
use warnings;
use Carp qw (cluck confess croak);
use Pod::Usage;
use Getopt::Long;
use File::Basename;

use lib '/NAS/arcad_data/Softs/tags';

use Modules::Config::Softwares;
use Modules::Files::Bam;

sub help
{ 
	my $msg = shift;
	
	print "-----------------------------------------------------------------------------\n"
	."-- ERROR MESSAGE --\n"."$msg\n".
	"-----------------------------------------------------------------------------\n"
	if( defined  $msg );
	pod2usage(0);
}

=head1 OPTIONS

=head2 required parameters

=over

=item B<-I|--input> (input file)

Specify your BAM file in input. Use this parameter as often as necessity.

=item B<-R|--reference> (reference file)

Reference used for mapping.

=item B<-knownSites> (already known variants)

A database of known polymorphic sites to skip over in the recalibration algorithm.
Use for BaseRecalibrator

=item B<-o|--out> (output file)

An existing directory in which to put all the results

=back

=head2 optional arguments

=over

=item B<-L | --intervals>

One or more genomic intervals over which to operate. 
Can be explicitly specified on the command line or in a file (including a rod file). 
Using this option one can instruct the GATK engine to traverse over only part of the genome. 
This argument can be specified multiple times. 
One may use samtools-style intervals either explicitly (e.g. -L chr1 or -L chr1:100-200) or listed in a file (e.g. -L myFile.intervals). 
Additionally, one may specify a rod file to traverse over the positions for which there is a record in the file (e.g. -L file.vcf). 
To specify the completely unmapped reads in the BAM file (i.e. those without a reference contig) use -L unmapped. 

=item B<-known>

Input VCF file with known indels. 
Any number of VCF files representing known SNPs and/or indels. 
Could be e.g. dbSNP and/or official 1000 Genomes indel calls. 
SNPs in these files will be ignored unless the --mismatchFraction argument is used. 
--known binds reference ordered data.
This argument supports ROD files of the following types: BCF2, VCF, VCF3

=cut

my $version="";

my  $bams = undef;
my(@opt_I, $opt_reference, $opt_o, @opt_known, @opt_knownSites, @opt_L, $opt_maxiter, $opt_queue, $fixMisencodedQuals);
my $min_variant_score = 10000;
my $keep_intermediate = 0;
#-----------------------------------------------------
# CREATION of output folders
#-----------------------------------------------------
my $realignment_dir   = 'real';
my $selection_dir     = 'select';
my $recalibration_dir = 'recal';
#-----------------------------------------------------

my %GATK_Input = (
	'-I'      => \@opt_I,
	'-R'      => \$opt_reference,
	'-L'      => \@opt_L,
);

my %GATK_RealignerTargetCreator = (
	'-o'      => sub{ mkdir "$opt_o/$realignment_dir/"; return "$opt_o/$realignment_dir/RTC.intervals";},
	'-known'  => \@opt_known,
	'-T'      => 'RealignerTargetCreator',
	%GATK_Input,
);

my %GATK_IndelRealigner = (
#	'-o' => sub{ my $tmp = $opt_o; $tmp =~ s/\.bam$//; return "$tmp.realigned.bam";},
	'-o' => sub{ return "$opt_o/$realignment_dir/realigned.bam";},
	'-known'  => \@opt_known,
	'-T'      => 'IndelRealigner',
	'-targetIntervals' => $GATK_RealignerTargetCreator{'-o'},
	%GATK_Input,
);

my %GATK_BaseRecalibrator = (
	'-o' => sub{ my $folder = "$opt_o/$recalibration_dir/"; mkdir $folder; return "${folder}BR.recal.grp";},
	'-knownSites' => \@opt_knownSites,
	'-T'          => "BaseRecalibrator",
#	'-plots'      => sub{ my $folder = "$opt_o/$recalibration_dir/"; return "${folder}BR.recal.plot.pdf";},
	'-I'          => $GATK_IndelRealigner{-o},
	'-R'          => \$opt_reference,
	'-L'          => \@opt_L,
#	%GATK_Input,
);

my %GATK_PrintReads = (
	'-o'    => sub{ my($file, $folder, $ext) = fileparse(&{$GATK_IndelRealigner{-o}}, qr/\.[^.]*/); my($fake, $folder2, $fake2) = fileparse(&{$GATK_BaseRecalibrator{-o}}, qr/\.[^.]*/); return "$folder2/${file}_recalibrated$ext";},
	'-BQSR' => $GATK_BaseRecalibrator{'-o'},
	'-T'    => 'PrintReads',
	'-I'    => $GATK_IndelRealigner{-o},
	'-R'    => \$opt_reference,
	'-L'    => \@opt_L,
#	%GATK_Input,
);

GetOptions(
	"help|?|h" => sub{ help(); },
	"input|I=s@"     	=> \@opt_I,
	"R|reference=s"  	=> \$opt_reference,
	"knownSites=s@"  	=> \@opt_knownSites,
	"o|out=s"        	=> \$opt_o,
	"known:s@"       	=> \@opt_known,
	"L|intervals:s"  	=> \@opt_L,
	"q|queue:s"      	=> \$opt_queue,
	"fixMisencodedQuals!"	=> \$fixMisencodedQuals
#	"Miter:i"        => \$opt_maxiter,
#	"keep_intermediate" => \$keep_intermediate,
);

help("No valid output found.\nPlease set the parameter '-o'") unless( $opt_o );
help("No valid input.\nPlease set the parameter '-i'") unless( @opt_I );
help("No valid reference. Please set the parameter '-r'") unless( $opt_reference && -e $opt_reference );

mkdir $opt_o if( defined $opt_o && $opt_o ne "" && !-d $opt_o );

my $o_unified = sub
{
	my($file, $folder, $ext) = fileparse($GATK_IndelRealigner{-o}, qr/\.[^.]*/); 
	$folder .= "/$selection_dir"; mkdir $folder;
	return "$folder/${file}_realigned.vcf";
};

#-----------------------------------------------------
# EXECUTABLES path
#-----------------------------------------------------
my $java_opts = 'Xmx5g'; # Pour les très gros fichiers, option inutile pour l'instant
my $JAVA_PATH = &$Softwares::JAVA_PATH or confess("$!");
my $GATK_DIR  = &$Softwares::GATK_DIRECTORY or confess("$!");
my $GATK_COMMAND = "$JAVA_PATH -jar $GATK_DIR/GenomeAnalysisTK.jar";
my $PICARD_TOOLS_DIRECTORY=&$Softwares::PICARD_TOOLS_DIRECTORY;
#-----------------------------------------------------

#-----------------------------------------------------
# INDEXING input bam file
#-----------------------------------------------------
my $o_bam = Bam->new( \@opt_I );
$o_bam->indexBam;
#-----------------------------------------------------

#--------------------------------------------------------------------------
# Création des fichiers : reference.fasta.fai et reference.dict s'ils n'existent pas 
#---------------------------------------------------------------------------

	if(! -e "$opt_reference.fai")
	{
		system("samtools faidx $opt_reference");
	}
	my $reference = $opt_reference;
	$reference =~ s/\.fasta//;
	$reference =~ s/\.fas//;
	$reference =~ s/\.fa//;
	if(! -e "$reference.dict")
	{
		my $command = "$JAVA_PATH -jar $PICARD_TOOLS_DIRECTORY/CreateSequenceDictionary.jar R=$opt_reference O=$reference.dict ";
		#print "\t$command\n";
		system($command);
	}
#----------------------------------------------------------------------------
if($fixMisencodedQuals)
{
	$GATK_RealignerTargetCreator{'-fixMisencodedQuals'}='';
	$GATK_IndelRealigner{'-fixMisencodedQuals'}='';
}
qsub_command( &make_command(\%GATK_RealignerTargetCreator, $GATK_COMMAND), 0 );
qsub_command( &make_command(\%GATK_IndelRealigner, $GATK_COMMAND), 0 );
Bam->new( [$GATK_IndelRealigner{-o}] )->indexBam;

#@opt_I = ( &{$GATK_IndelRealigner{-o}} );
if( @opt_knownSites )
{
	qsub_command( &make_command(\%GATK_BaseRecalibrator, $GATK_COMMAND) );
	qsub_command( &make_command(\%GATK_PrintReads, $GATK_COMMAND) );
	Bam->new( [$GATK_PrintReads{-o}] )->indexBam;
}
else
{
	#Creer un fichier de recalibration
#	my %unified = (-I => \@opt_I, -R => $GATK_Input{-R}, -o => "$opt_o.all.vcf");
	my %unified = (
		-I => $GATK_PrintReads{-I}, 
		-R => $GATK_Input{-R}, 
		-o => sub
		{
			my($file, $folder, $ext) = fileparse(&{$GATK_IndelRealigner{-o}}, qr/\.[^.]*/); 
			$folder = "$opt_o/$selection_dir"; mkdir $folder;
			return "$folder/${file}_realigned.vcf";
		}
	);
	$unified{-L} = $GATK_Input{-L} if( exists $GATK_Input{-L} );
	my $genotyper_call_script = &$Softwares::ARCAD_SCRIPT_PATH('7_genotyper_call.pl');
	#qsub_command( &make_command(\%unified, "perl /home/homa/arcad/Softs/trunk/7_genotyper_call.pl --filtration --nophasing --standard_min_confidence_threshold_for_calling $min_variant_score --standard_min_confidence_threshold_for_emitting $min_variant_score"), 0);
	qsub_command( &make_command(\%unified, "perl $genotyper_call_script --filtration --nophasing --standard_min_confidence_threshold_for_calling $min_variant_score --standard_min_confidence_threshold_for_emitting $min_variant_score"), 0);
	
	# Filter PASS from vcf
#	my $selected_vcf = "$opt_o.selected.tmp.vcf"; #$min_variant_score
#	`awk '/^#/{print \$0 > \"$selected_vcf\"} {if(\$6 >= $min_variant_score && \$7==\"PASS\"){print \$0 >> \"$selected_vcf\"}}' "$opt_o.all.vcf_filtered.vcf"`;
	
	my %GATK_SelectVariants = (
		'-T'        => 'SelectVariants',
		'-V'        => &{$unified{-o}} . "_filtered.vcf",
		'-selectType' => 'SNP',
		'-select'   => '"vc.isNotFiltered()"',
		'-R'        => \$opt_reference,
		'-L'        => \@opt_L,
		'-o'        => sub
		{
			my($file, undef, undef) = fileparse(&{$GATK_IndelRealigner{-o}}, qr/\.[^.]*/); 
			my $folder = "$opt_o/$selection_dir"; #mkdir $folder;
			return "$folder/${file}_selected.vcf";
		}
	);
	qsub_command( &make_command(\%GATK_SelectVariants, $GATK_COMMAND), 0);
	
	push(@opt_knownSites, &{$GATK_SelectVariants{-o}});
	qsub_command( &make_command(\%GATK_BaseRecalibrator, $GATK_COMMAND) );
	qsub_command( &make_command(\%GATK_PrintReads, $GATK_COMMAND) );
	Bam->new( [$GATK_PrintReads{-o}] )->indexBam;
}

if( !$keep_intermediate )
{
	#&erase_intermediate_files;
}

sub erase_intermediate_files
{
	unlink $GATK_RealignerTargetCreator{-o} if( -e $GATK_RealignerTargetCreator{-o} );
	unlink $GATK_IndelRealigner{-o}         if( -e $GATK_IndelRealigner{-o} );
	unlink $GATK_BaseRecalibrator{-o}       if( -e $GATK_BaseRecalibrator{-o} );
	
	if( @opt_knownSites )
	{
		unlink "$opt_o.all.vcf_filtered.vcf" if( -e "$opt_o.all.vcf_filtered.vcf" );
		unlink "$opt_o.all.vcf"              if( -e "$opt_o.all.vcf" );
		unlink "$opt_o.selected.tmp.vcf"     if( -e "$opt_o.selected.tmp.vcf" );
		unlink "$opt_o.selected.vcf"         if( -e "$opt_o.selected.vcf" );
	}
}

sub make_command
{
	my $rh_command = shift;
	my $prefix     = shift;
	my $suffix     = shift;
	
	my $command = ($prefix)? "$prefix " : '';
	foreach my $key ( keys %$rh_command )
	{
		if( ref $$rh_command{$key} eq 'ARRAY' )
		{
			my $arr = $$rh_command{$key};
			foreach my $val ( @$arr )
			{
				$command .= "$key $val ";
			}
		}
		elsif( ref $$rh_command{$key} eq 'SCALAR' )
		{
			$command .= "$key ${$$rh_command{$key}} " if( defined ${$$rh_command{$key}} );
		}
		elsif( ref $$rh_command{$key} eq 'CODE' )
		{
			my $tmp = &{$$rh_command{$key}};
			$command .= "$key $tmp ";
		}
		else
		{
			$command .= "$key $$rh_command{$key} ";
		}
	}
	return ($suffix)? "$command $suffix" : $command;
}

sub qsub_command
{
	my $command = shift;
	my $debug = shift;
	
	if( $ENV{JOB_ID} )
	{
		if($debug){
			print "$command\n";
		}
		else
		{
			$? = system( $command );
			help("Job '$ENV{JOB_NAME}' with job id '$ENV{JOB_ID}' failed") if( $? != 0 );
		}
	}
	else
	{
		my $rh_qsub = &$Softwares::RH_QSUB; 
		$$rh_qsub{q} = $opt_queue if( $opt_queue );
		my $qsub = &Softwares::make_qsub_command( $rh_qsub );
		if($debug){
			print "$qsub '$command'\n";
		}
		else{
			$? = system( "$qsub '$command'" );
			help("Unable to run command \n'$qsub $command'\n") if( $? != 0 );
		}
	}
}

=head1 DEPENDENCIES

=head1 INCOMPATIBILITIES

<Fully compatble with all version of perl :)>

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 1

=item Felix Homa (Cirad) felixhoma@yahoo.fr

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

=cut
