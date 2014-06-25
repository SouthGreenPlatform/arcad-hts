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
=head1 NAME

GenotypeConcordance.pl - Compares pairs of VCF files 
and each vcf file with known SNPs(VCF file)if available

=head1 DESCRIPTION

Compares pairs of VCF files and each vcf file with known SNPs(VCF file) if available

Makes a statistical summary over a NGS/GBS data analysis in a file

Generates a graph of missing data for each vcf file, a report of the comparison and summary 

file (total number of SNPs and maximum value of missing data per individual and per SNPs, in each file)


=head1 SYNOPSIS / USAGE

GenotypeConcordance.pl -R reference -directory directory -comp Known SNPs(vcf file) -output prefix of the output folder

GenotypeConcordance.pl -R reference -d directory -c Known SNPs(vcf file) -o prefix_of_the_output_older

=cut

use strict;
use Carp qw (cluck confess croak);
#use warnings;
use Getopt::Long;
use Pod::Usage;
use Fatal qw/:void open close/;
use Cwd qw(cwd abs_path);
use Tie::File;
use List::Util qw(max);
use Statistics::R;

use lib '/NAS/arcad_data/Softs/tags';
use Modules::Config::Softwares;
use Modules::Files::Files;



=head1 OPTIONS

=head2 required parameters

=over 5

=item B<[-R]> ([reference file]):

Fasta reference file

=item B<[-d]> ([a directory]): 

Directory containing VCF files to evaluate

=item B<[-comp]> ([vcf file]):

Variants and genotypes to compare against (vcf file)

=item B<[-o]> ([prefix output directory]):

Folder where statistics, for each comparaison, will be written



=back

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

pod2usage(0) unless (@ARGV);


#options processing
my ($man, $help, $version, $ref, $eval,$directory, $comp, $output,$vcf_files, @opt_I, $out);
my $pattern = '*.vcf';
my $subdirectory = 0;

my $QUEUE = 'arcad.q';

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"R|reference=s" => \$ref,
	"directory|d=s"    => \$directory,
	"comp|c=s"    => \$comp,
	"output|o=s" => \$output,
	"subdirectory!" => \$subdirectory,
	"pattern|p=s" => \$pattern
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

help "No VCF file(s) found" unless( $directory );

if( -d $directory )
{
	$vcf_files = Files->getFiles( $directory, $pattern, $subdirectory );
	help "No VCF file(s) found" if( 0 >= scalar @$vcf_files );
}
else{
	help "$! $directory";
}


#-----------------------------------------------------
# EXECUTABLES path
#-----------------------------------------------------
my $java_opts = 'Xmx5g'; 
my $JAVA_PATH = &$Softwares::JAVA_PATH or confess("$!");
my $GATK_DIR  = &$Softwares::GATK_DIRECTORY or confess("$!");
my $GATK_COMMAND = "$JAVA_PATH -jar $GATK_DIR/GenomeAnalysisTK.jar";
#-----------------------------------------------------

#create the output folder in the current directory 
my $tmp = $output; 

my ($nbr_SNP, $vcf, $out1, %hash, $taille, @tab, @tableau, @tableau2, @b, @a, $report, $jpg, $path, $files_compare, @vcf_filess);
mkdir($tmp);

$report =$output."/".$output."_report.txt";
open (OUT, ">$report") or die "cannot open the file $report\n";

foreach $vcf ( @$vcf_files )
{
	my $NN = 0;
	@b = split('/', $vcf);
	$out1 = $output."/".$b[$#b]."_NA.txt";
	open (OUT1, ">$out1") or die "cannot open the file $out1\n";
	open (F, $vcf) or die "cannot open the file $vcf\n";
	$nbr_SNP = `grep -c '^[^#]' $vcf`;
	print OUT "Nombre de SNPs ($vcf) : $nbr_SNP";
	push (@opt_I,$vcf);
	# generate a matrix with presence "1" absence "0"of data,to generate aplot of missing data with R
	while(<F>){
		chomp($_);
		@tab = split(/\t/, $_);
		if($_ =~ /^#/)
		{
			if($_ =~ /^#CHROM/)
			{
				print OUT1 "SNP\t";
				for(my $i = 9;$i<@tab;$i++){
					if($i== $#tab){
						print OUT1 "$tab[$i]";
					}
					else{
						print OUT1 "$tab[$i]\t";
					}
				}
				print OUT1 "\n";
			}
		}
		else{
			print OUT1 "$tab[0]_$tab[1]\t";
			for(my $i = 9;$i<@tab;$i++){
				if ($tab[$i] =~ m/^(\.)\/(\.)/){
					if($i== $#tab){
						print OUT1 "0";
					}
					else{
						print OUT1 "0\t";
					}
					$NN++;
					if ( defined $hash{$i} ){
						$hash{$i}++;
					}
					else{
						$hash{$i} = 1;
					}
				}
				else{
					if($i== $#tab){
						print OUT1 "1";
					}
					else{
						print OUT1 "1\t";
					}
				}
			}
			print OUT1 "\n";
			push(@tableau, $NN);
			$NN = 0;
		}
	}

for( my $j=9;$j<@tab;$j++){
	
    if ( defined $hash{$j} ){
        #print OUT1 "$hash{$j}\t";
        push(@tableau2, $hash{$j});
    }
    else{
        #print OUT1 "0\t";
    }
}


my $INDIV = @tab - 9;
my $NA_indiv = max(@tableau) / $INDIV * 100;
my $NA_snp = max(@tableau2) / $nbr_SNP * 100;
$nbr_SNP =~ s/\n//;
print OUT "NA/INDIV :".max(@tableau)."/".$INDIV." soit : $NA_indiv % \n";
print OUT "NA/SNP : ".max(@tableau2)."/".$nbr_SNP." soit : $NA_snp % \n";

@tableau = 0;
@tableau2 = 0;
$NN = 0;
%hash = ();

#-----------------------------------------------------
# R : Missing data graph
#-----------------------------------------------------
			my $r_script= "GenotypeConcordance2.R";
			my $path="./$r_script";
			$jpg = $output."/".$b[$#b].".png";
			$jpg =~ s/\.vcf//;
			my $execute = `Rscript $path -n $out1 -jpg $jpg `;
			print $? if $?;
			print $execute,"\n";
#-----------------------------------------------------


if($comp ne ''){
#-----------------------------------------------------
# Genotype Concordance 
#-----------------------------------------------------

			@a = split('/', $comp);
			$out = $output."/".$b[$#b]."_Vs_".$a[$#a]."_GenotypeConcordanceReport.report";
			$report = $output."/".$b[$#b]."_Vs_".$a[$#a]."_Report.txt";
			if(-e ! $report){
				open (R, ">$report") or die "cannot open the file $report\n";
			}
			my $genocom1 = "$GATK_COMMAND -T GenotypeConcordance -R $ref -o $out -eval $vcf -comp $comp --ignoreFilters";
			print R "file1 = Eval = $vcf \n";
			print R "file2 = Comp = $comp \n";
			print " $genocom1\n\n";
			system ($genocom1 );
			sleep(10);
#-------------------------------------------
# R : Read GATK report
#-------------------------------------------
			my $r_script= "GenotypeConcordance.R";
			my $path="./$r_script";
			$jpg = $output."/".$b[$#b].".png";
			$jpg =~ s/\.vcf//;
			my $execute = `Rscript $path -i $out -r $report `;
			print $? if $?;
			print $execute,"\n";
			sleep(10);
#-------------------------------------------
#-----------------------------------------------------
}
else{
#-----------------------------------------------------
# Genotype Concordance 
#-----------------------------------------------------

	for(my $i=0; $i<@$vcf_files; $i++){
		if($vcf ne @$vcf_files[$i]){
			$files_compare = "-eval $vcf -comp @$vcf_files[$i]";
			@a = split('/', @$vcf_files[$i]);
			$out = $output."/".$b[$#b]."_Vs_".$a[$#a]."_GenotypeConcordanceReport.report";
			$report = $output."/".$b[$#b]."_Vs_".$a[$#a]."_Report.txt";
			if(!-e $report){
				open (R, ">$report") or die "cannot open the file $report\n";
			}
			my $genocom1 = "$GATK_COMMAND -T GenotypeConcordance -R $ref -o $out $files_compare --ignoreFilters ";
			
			print R "file1 = Eval = $vcf \n";
			print R "file2 = Comp = @$vcf_files[$i] \n";
			
			print " $genocom1\n\n";
			system ($genocom1 );
			
#-------------------------------------------
# R : Read GATK report
#-------------------------------------------
			my $r_script= "GenotypeConcordance.R";
			my $path="./$r_script";
			$jpg = $output."/".$b[$#b].".png";
			$jpg =~ s/\.vcf//;
			my $execute = `Rscript $path -i $out -r $report `;
			print $? if $?;
			print $execute,"\n";
#-------------------------------------------
		}
	}
#-----------------------------------------------------
}
system("rm $out");
}


system("rm $out");


=head1 DEPENDENCIES

GATK

=head1 INCOMPATIBILITIES

<Fully compatble with all version of perl :)>

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 1

=item Hajar CHOUIKI (INRA) hajar.chouiki@supagro.inra.fr

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
