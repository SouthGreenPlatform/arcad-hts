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

=head1 genotyper_call.pl

genotyper_call.pl - Use the GATK software to make variant call in NGS datat processing

=head1 DESCRIPTION

This script use UnifiedGenotyper module of the GATK software to make variant call.
Parameter are the same than GATK 'SHORT PARAMETER'.

=head1 SYNOPSIS / USAGE

[qsub -b y -cwd -q arcad.q -N variantsCaller] genotyper_call.pl -I bam_file.bam -o out.vcf -R reference.fasta

[qsub -b y -cwd -q arcad.q -N variantsCaller] genotyper_call.pl -I bam1.bam -I bam2.bam -I bam3.bam -o out.vcf -R reference.fasta

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

# GATK VARIANTS FILTRATION
our $opt_dp_filter = 20;
our $opt_clusterSize = 3;
our $opt_clusterWindowSize = 10;
our @opt_filterExpression = ('"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"', '"QD < 1.5"', '"DP < 20"');
our @opt_filterName = ("HARD_TO_VALIDATE", "QD_FILTER",  "DP_FILTER");
our @opt_genotypeFilterExpression;
our @opt_genotypeFilterName;
our $opt_invalidatePreviousFilters = 'false';
our $opt_mask;
our $opt_maskExtension = 0;
our $opt_maskName = 'Mask';
our $opt_missingValuesInExpressionsShouldEvaluateAsFailing = 'false';

# Input/output files
our $opt_variant;
our $opt_out;
our $opt_o;
our @opt_I = ();
our $opt_R;

my @opt_intervals;

# GATK READ BACKED PHASING
our $opt_cacheWindowSize = 20000;
our $opt_enableMergePhasedSegregatingPolymorphismsToMNP;
our $opt_maxGenomicDistanceForMNP = 1;
our $opt_maxPhaseSites = 10;
our $opt_min_base_quality_score = 17;
our $opt_min_mapping_quality_score = 20;
our $opt_phaseQualityThresh = 10.0;
our $opt_respectPhaseInInput;
our $opt_sampleToPhase;

=item

our %GATK_IN_OUT = (
	variant  => \$opt_variant,
	V        => \$opt_variant,
	out      => \$opt_out,
	o        => \$opt_out,
);

=cut

my %GATK_STANDARD = (
	'intervals' => \@opt_intervals,
);

#my $opt_dp_filter = undef;
my %h_GATK_VARIANTSFILTRATION_LONG = (
	clusterSize                                       => \$opt_clusterSize,
	clusterWindowSize                                 => \$opt_clusterWindowSize,
	filterExpression                                  => \@opt_filterExpression,
	filterName                                        => \@opt_filterName,
	genotypeFilterExpression                          => \@opt_genotypeFilterExpression,
	genotypeFilterName                                => \@opt_genotypeFilterName,
	invalidatePreviousFilters                         => \$opt_invalidatePreviousFilters,
	mask                                              => \$opt_mask,
	maskExtension                                     => \$opt_maskExtension,
	maskName                                          => \$opt_maskName,
	missingValuesInExpressionsShouldEvaluateAsFailing => \$opt_missingValuesInExpressionsShouldEvaluateAsFailing,
#	%GATK_STANDARD
);

our %h_GATK_READBACKEDPHASING_LONG = (
	cacheWindowSize                                       => \$opt_cacheWindowSize,
	enableMergePhasedSegregatingPolymorphismsToMNP        => \$opt_enableMergePhasedSegregatingPolymorphismsToMNP,
	maxGenomicDistanceForMNP                              => \$opt_maxGenomicDistanceForMNP,
	maxPhaseSites                                         => \$opt_maxPhaseSites,
	min_base_quality_score                                => \$opt_min_base_quality_score,
	min_mapping_quality_score                             => \$opt_min_mapping_quality_score,
	phaseQualityThresh                                    => \$opt_phaseQualityThresh,
	respectPhaseInInput                                   => \$opt_respectPhaseInInput,
	sampleToPhase                                         => \$opt_sampleToPhase,
	%GATK_STANDARD
);

my $opt_alleles;
#my $opt_annotateNDA = 'false';
my $opt_annotateNDA;
my @opt_annotation;
my @opt_comp;
#my $opt_computeSLOD = 'false';
my $opt_computeSLOD;
#my $opt_contamination_fraction_to_filter = 0.05;
my $opt_contamination_fraction_to_filter;
my @opt_excludeAnnotation;
my $opt_genotype_likelihoods_model = 'BOTH';
my $opt_genotyping_mode = 'DISCOVERY';
#my @opt_group = ('Standard');
my @opt_group;
my $opt_heterozygosity = 0.001;
#my $opt_ignoreLaneInfo = 'false';
my $opt_ignoreLaneInfo;
my $opt_output_mode = 'EMIT_VARIANTS_ONLY';
my $opt_dfrac = 1;
my $opt_standard_min_confidence_threshold_for_calling = 50.0;
my $opt_standard_min_confidence_threshold_for_emitting = 30.0;
my $opt_downsampling_type = 'ALL_READS';

my %h_GATK_UNIFIEDGENOTYPER_LONG = (
	'alleles' => \$opt_alleles,
	'annotateNDA' => \$opt_annotateNDA,
	'annotation' => \@opt_annotation,
	'comp' => \@opt_comp,
	'computeSLOD' => \$opt_computeSLOD ,
	'contamination_fraction_to_filter' => \$opt_contamination_fraction_to_filter,
	'excludeAnnotation' => \@opt_excludeAnnotation,
	'genotype_likelihoods_model' => \$opt_genotype_likelihoods_model,
	'genotyping_mode' => \$opt_genotyping_mode,
	'group' => \@opt_group,
	'heterozygosity' => \$opt_heterozygosity,
	'ignoreLaneInfo' => \$opt_ignoreLaneInfo,
	'output_mode' => \$opt_output_mode,
	'downsample_to_fraction' => \$opt_dfrac,
	'standard_min_confidence_threshold_for_calling' => \$opt_standard_min_confidence_threshold_for_calling,
	'standard_min_confidence_threshold_for_emitting' => \$opt_standard_min_confidence_threshold_for_emitting,
	'downsampling_type' => \$opt_downsampling_type,
	%GATK_STANDARD
);


sub help
{ 
	my $msg = shift;
	
	print "-----------------------------------------------------------------------------\n"
	."-- ERROR MESSAGE --\n"."$msg\n".
	"-----------------------------------------------------------------------------\n"
	if( defined  $msg );
	pod2usage(0);
}

help 'No parameter detected, use -h option to see help' unless (@ARGV);

=head1 OPTIONS

=head2 required parameters

=over

=item B<-I> (input file)

Specify your BAM file in input. Use this parameter as often as necessity.

=item B<-R> (Reference file)

Reference used for mapping

=item B<-o> (output file)

File in which first vcf will be write, it will be reused as pattern for 
filtration and phasing output file name

=item B<--filtration> (boolean, default is false)

--filtration, if you want to apply the variantFiltration on data
--nofiltration, if you don't want to apply the variantFiltration on data

=item B<--phasing> (boolean, default is false)

Specify if you want to make the phasing step in the analysis by specifying --phasing, 
or --nophasing if you don't want to phase data


=back

=head2 optional arguments

=over

=item B<all gatk unifiedGenotyper options>

Check here for all possibility http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

=back

=cut

my $version="";

my($phasing, $filtration);

my $QUEUE = 'arcad.q';
my $JOBNAME = 'unifiedgenotyper';
my $LOGFILE = "variantsCall.o$$";
my $QSUB = "qsub -b y -cwd -q <<Q>> -N <<JOBNAME>> -sync y -o $LOGFILE ";

GetOptions(
#	initial parameters
	"phasing!" => \$phasing,
	"filtration!"  => \$filtration,
	"help|?|h" => sub{ help(); },
#	Input
	'I=s@',
	'o|out=s',
	'R|reference=s' => \$opt_R,
	
	'L|intervals:s@' => \@opt_intervals,
	
#	ReadBackPhasing
	'cacheWindow|cacheWindowSize:i',
	'enableMergeToMNP|enableMergePhasedSegregatingPolymorphismsToMNP:s',#bool
	'maxDistMNP|maxGenomicDistanceForMNP:i',
	'maxSites|maxPhaseSites:i',
	'mbq|min_base_quality_score:i',
	'mmq|min_mapping_quality_score:i',
	'phaseThresh|phaseQualityThresh:i',
	'respectPhaseInInput:s',#bool
	'sampleToPhase:s',
	
#	Variant Filtration
	'clusterSize|cluster:i' => \$opt_clusterSize,
	'window|clusterWindowSize:i'  => \$opt_clusterWindowSize,
	'filter|filterExpression:s@' => \@opt_filterExpression,
	'depth:i'                    => \$opt_dp_filter,
	'filterName:s@',
	'G_filter|genotypeFilterExpression:s@',
	'G_filterName|genotypeFilterName:s@',
	'invalidatePreviousFilters:s',#bool
	'mask:s',
	'maskExtend|maskExtension:i',
	'maskName:s',
	'missingValuesInExpressionsShouldEvaluateAsFailing:s',#bool
	
#	Unified Genotyper
	'alleles:s' => \$opt_alleles,
	'annotateNDA:s' => \$opt_annotateNDA,
	'annotation:s@' => \@opt_annotation,
	'comp:s@' => \@opt_comp,
	'computeSLOD:s' => \$opt_computeSLOD ,
	'contamination_fraction_to_filter:f' => \$opt_contamination_fraction_to_filter,
	'excludeAnnotation:s@' => \@opt_excludeAnnotation,
	'genotype_likelihoods_model:s' => \$opt_genotype_likelihoods_model,
	'genotyping_mode:s' => \$opt_genotyping_mode,
	'group:s@' => \@opt_group,
	'heterozygosity:f' => \$opt_heterozygosity,
	'ignoreLaneInfo:f' => \$opt_ignoreLaneInfo,
	'output_mode:s' => \$opt_output_mode,
	'downsample_to_fraction:f' => \$opt_dfrac,
	'standard_min_confidence_threshold_for_calling:f' => \$opt_standard_min_confidence_threshold_for_calling,
	'standard_min_confidence_threshold_for_emitting:f' => \$opt_standard_min_confidence_threshold_for_emitting,
	'downsampling_type:s' => \$opt_downsampling_type,
);

help 'No valid BAM file(s)' unless ( 0 < scalar @opt_I );
help 'No valid REFERENCE file *'.$opt_R.'*' unless ( $opt_R && -e $opt_R );
help 'Please set OUTPUT parameter (-o output_file)' unless ( $opt_o );
$opt_filterExpression[2] = '"DP < ' . $opt_dp_filter . '"' if( defined $opt_dp_filter );


my $filter_out = "${opt_o}_filtered.vcf";
my $phase_out  = "${opt_o}_phased.vcf";

#-----------------------------------------------------
# EXECUTABLES path
#-----------------------------------------------------
my $java_opts = 'Xmx4g'; # Pour les tr\E8s gros fichiers, option inutile pour l'instant
my $JAVA_PATH = &$Softwares::JAVA_PATH or confess("$!");
my $GATK_DIR  = &$Softwares::GATK_DIRECTORY or confess("$!");
my $GATK_COMMAND = "$JAVA_PATH $java_opts -jar $GATK_DIR/GenomeAnalysisTK.jar";
my $SAMTOOLS_PATH = &$Softwares::SAMTOOLS_PATH or confess("$!");
#-----------------------------------------------------

sub makeCommands
{
	my %h_h = @_;
	
	my $command_ext = '';
	foreach my $key ( keys %h_h )
	{
		if( ref $h_h{$key} eq 'ARRAY' )
		{
			next if( 0 >= @{$h_h{$key}} );
			$command_ext .= "--$key $_ " foreach ( @{$h_h{$key}} );
		}
		else
		{
			next if( !defined ${$h_h{$key}} );
			next if( ${$h_h{$key}} =~ m/^\s*false\s*$|^f$/i );
			if( ${$h_h{$key}} =~ m/^\s*true\s*$|^t$/i )
			{
				$command_ext .= "--$key ";
			}
			else{
				$command_ext .= "--$key ${$h_h{$key}} ";
			}
		}
	}
	return $command_ext;
}

#-----------------------------------------------------
# INDEXING bam file
#-----------------------------------------------------
my $o_bam = Bam->new( \@opt_I );
$o_bam->indexBam;
#-----------------------------------------------------

#-----------------------------------------------------
# UNIFIED GENOTYPER : VARIANTS CALL
#-----------------------------------------------------
my $input_command = '-I ' . join(' -I ', @opt_I);
my $genocom = "$GATK_COMMAND -T UnifiedGenotyper " . &makeCommands( %h_GATK_UNIFIEDGENOTYPER_LONG ) . " $input_command -o $opt_o -R $opt_R";
print "$genocom\n\n";
if( $ENV{JOB_ID} )
{
	$? = system( $genocom );
	help("Job '$ENV{JOB_NAME}' with job id '$ENV{JOB_ID}' failed") if( $? != 0 );
}
else
{
	my $qsub = $QSUB;
	$qsub =~ s/\<\<Q\>\>/$QUEUE/;
	$qsub =~ s/\<\<JOBNAME\>\>/unified_genotyper/;
	$? = system( "$qsub $genocom" );
	help("Job 'unified_genotyper' failed") if( $? != 0 );
}
#-----------------------------------------------------

#-----------------------------------------------------
# VARIANTS FILTRATION
#-----------------------------------------------------
if( $filtration )
{
	my $genocom = "$GATK_COMMAND -T VariantFiltration " . &makeCommands( %h_GATK_VARIANTSFILTRATION_LONG ) . " -V $opt_o -R $opt_R --out $filter_out";
	print "$genocom\n\n";
	if( $ENV{JOB_ID} )
	{
		$? = system( $genocom );
		help("Job '$ENV{JOB_NAME}' with job id '$ENV{JOB_ID} failed") if( $? != 0 );
	}
	else
	{
		my $qsub = $QSUB;
		$qsub =~ s/\<\<Q\>\>/$QUEUE/;
		$qsub =~ s/\<\<JOBNAME\>\>/variant_filtration/;
		$? = system( "$qsub '$genocom'" );
		help("Job 'variant_filtration' failed") if( $? != 0 );
	}
}
#-----------------------------------------------------

#-----------------------------------------------------
# READ BACKED PHASING
#-----------------------------------------------------
if( $phasing )
{
	my $phase_input = ($filtration)? $filter_out:$opt_o;
	my $genocom = "$GATK_COMMAND -T ReadBackedPhasing " . &makeCommands( %h_GATK_READBACKEDPHASING_LONG ) . " $input_command -V $phase_input -R $opt_R --out $phase_out";
	print "$genocom\n\n";
	if( $ENV{JOB_ID} )
	{
		$? = system( $genocom );
		help("Job '$ENV{JOB_NAME}' with job id '$ENV{JOB_ID} failed") if( $? != 0 );
	}
	else
	{
		my $qsub = $QSUB;
		$qsub =~ s/\<\<Q\>\>/$QUEUE/;
		$qsub =~ s/\<\<JOBNAME\>\>/phasing/;
		$? = system( "$qsub $genocom" );
		help("Job 'phasing' failed") if( $? != 0 );
	}
}
#-----------------------------------------------------

=head1 DEPENDENCIES

GATK 2.0

=head1 INCOMPATIBILITIES

<Fully compatble with all version of perl :)>

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 1

=item Felix HOMA (IRD) felix.homa-at-supagro.inra.fr

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
