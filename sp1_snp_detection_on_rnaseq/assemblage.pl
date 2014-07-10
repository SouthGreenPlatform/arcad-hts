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

assemblage.pl - Perform sequences assembly

=head1 DESCRIPTION

Perform sequences assembly using one or several assemblers

Currently, you can use abyss-pe and cap3.

=head1 SYNOPSIS

assemblage.pl 
[-i|-input input_folder_directory] 
[-o|-output output_folder] 
[-n|-name analysis_name] 
[-f|-forward forward_pattern] 
[-r|-reverse reverse_pattern] 
[-abyss abyss_options] 
[-cap3 cap3_options] 

for example: I want to run assemblage on paired end reads (forward, reverse) conatained on folder 
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/

LIST OF FILES IN THE FOLDER E<10>
$ ls /NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/*.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC10_TAGCTT_R1.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC10_TAGCTT_R2.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC10_TAGCTT_Single.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC1_ATCACG_R1.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC1_ATCACG_R2.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC1_ATCACG_Single.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC2_CGATGT_R1.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC2_CGATGT_R2.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC2_CGATGT_Single.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC3_TTAGGC_R1.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC3_TTAGGC_R2.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC3_TTAGGC_Single.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC4_TGACCA_R1.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC4_TGACCA_R2.fastq E<10>
/NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/Cleaned_FC4_TGACCA_Single.fastq
...

For Paired end reads
qsub -q bioinfo.q -N fonio_assembly -b y 
"perl /NAS/arcad_data/scripts/assemblage.pl -n fonio -i /NAS/arcad_data/sp1_final/fonio_millet/cleaned_data/ -o ~/work/assembly -f _R1 -r _R2 --abyss --cap3"

the Global pattern _R1 is specific of forward reads while _R2 is specific of reverse reads.

=cut

use strict;
use warnings;
use Carp qw (cluck confess croak);
use Readonly;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Getopt::Long;
use File::Path;
use Cwd;

use lib '/NAS/arcad_data/Softs/tags';
use Modules::Files::Files;
use Modules::Config::Softwares;

sub help
{ 
	my $msg = shift;
	
	#print "-----------------------------------------------------------------------------\n"
	#."-- ERROR MESSAGE --\n"."$msg\n".
	#"-----------------------------------------------------------------------------\n"
	#if( defined  $msg );
	
	pod2usage(
		-msg     => $msg,
		-exitval => 1,  
		-verbose => 2,
	);
}

help "Not enough arguments" unless (@ARGV);

=head1 OPTIONS

=head2 required parameters

=over 11

=item [-i|--input] <INPUT_FILE> or <INPUT_DIRECTORY>

Sequence file or directory containing sequences files. This parameter can be repetead.

=back

=head2 optional arguments

=over 11

=item B<[-o|--output]> <OUTPUT_DIRECTORY>

Save the results of analysis in the output directory. The default is your home directory.

=item B<[-n|--name]> <NAME>

Analysis name. The results of the current analysis will be stored in /OUTPUT_DIRECTORY/NAME. The default is 'ASSEMBLY'.

=item B<[-f|--forward]> <PERL_REGEXP_PATTERN>

PERL Pattern that will be used to retrieve all the file which contains forward straud. The default is 'forward'.

=item B<[-r|--reverse]> <PERL_REGEXP_PATTERN>

PERL Pattern that will be used to retrieve all the file which contains reverse straud. The default is 'reverse'.

=item B<[--cap3]> <' OPTIONS'>

Run cap3 assembler with specified options if any. Options must be surrounded by quotes and begin by a space. Exemple : assemblage.pl --cap3 ' -k 60'. By default, cap3 will be run without option.

=item B<[--abyss]> <' OPTIONS'>

Run abyss-pe assembler with specified options if any. Options must be surrounded by quotes and begin by a space. Exemple : assemblage.pl --abyss ' -k 60'. By default, cap3 will be run with the following options : name=assemblage_abyss OVERLAP_OPTIONS=--no-scaffold SIMPLEGRAPH_OPTIONS=--no-scaffold E=0 n=10 v=-v k=60.

=item B<[-h|-help|-?]>

Prints a brief help message and exits.

=item B<[-man]>

Prints detailled help and script documentation.

=back

=cut

my $number = "1"; # Numero d'etape

my $man;
my $help;
my $name = "assembly";
my @input = ();
my $output = $ENV{'HOME'};
my $forward_pattern = "forward";
my $reverse_pattern = "reverse";
my @analysis = ();

my $ABYSS_PE_PATH = &$Softwares::ABYSS_PE_PATH or confess("$!");
my $CAP3_PATH = &$Softwares::CAP3_PATH or confess("$!");

GetOptions(

			"help|?|h"			=> \$help,
			'man'				=> \$man,
			"n|name=s"			=> \$name,	
			"i|input=s" 		=> \@input,
			"o|output=s" 		=> \$output,
			"f|forward=s" 		=> \$forward_pattern,
			"r|reverse=s" 		=> \$reverse_pattern,
			"cap3:s" 		=> sub { register_analysis(\&cap3, $CAP3_PATH, $_[1]) },
			"abyss:s"	 	=> sub { register_analysis(\&abyss_pe, $ABYSS_PE_PATH, $_[1]) }

) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

pod2usage(0) if(!@input);

# List input files

if (!@input){

	@input = cwd();

}

my @file_list=();

foreach my $fichier (@input){

	if (-e $fichier){

		if (-d $fichier){

			push(@file_list, ls($fichier));


		}

		else{

			push(@file_list, $fichier);

		}

	}

}



# Running analysis



print "nb etapes = ".scalar(@analysis)."\n";

if (scalar(@analysis)){

print "Running\n";

	foreach my $step (@analysis){

		my $output_path = $output.'/'.$name.'/STEP'.$number.'_';
		@file_list = &{$step->[0]}($step->[1], \@file_list, $output_path, $step->[2]);
		$number++;

	}

}




sub register_analysis{

	my ($function_reference, $path, $parameters)=@_;
	push (@analysis, [$function_reference, $path, $parameters]);
}


sub cap3{

	my ($exec, $inputs, $output, $parameters) = @_;

	$output .= 'CAP3';

	mkpath($output);

	my $input=merge_files($output.'/input_cap3.fasta', @$inputs);

	my $commandLine = "$exec " . $input." $parameters ".' >'.$output.'/summary';

	print $commandLine."\n";	

	system($commandLine);


	prefix_fasta_sequence($output.'/input_cap3.fasta.cap.singlets', 'singlet ');

	return ls ($output, '.+\.(contigs|singlets)$') ;


}

sub abyss_pe{

	my ($exec, $inputs, $output, $parameters) = @_;

	if (! defined($parameters) || $parameters eq "")
	{
		$parameters=" name=assemblage_abyss OVERLAP_OPTIONS=--no-scaffold SIMPLEGRAPH_OPTIONS=--no-scaffold E=0 n=10 v=-v k=60"
	};

	$output .= 'ABYSS';

	mkpath($output);
	my @forward = grep {/$forward_pattern/} @$inputs;
	my @reverse = grep {/$reverse_pattern/} @$inputs;


	my $fd=merge_files($output.'/abyss_forward.fastq', @forward);
	my $rv=merge_files($output.'/abyss_reverse.fastq', @reverse);
	my $commandLine = "$exec $parameters" . ' in="'.$fd.' '.$rv.'" -C '.$output.' > '.$output.'/summary';

	print "$commandLine\n";
 
	system($commandLine);
	
	
	return ls ($output, '.+-contigs\.fa$');


}

sub merge_files{

	my ($name, @files_list) = @_;
	print 'cat '.join(" ",@files_list).' > '.$name."\n";
	system ('cat '.join(" ",@files_list).' > '.$name);
	
	return $name;

}


sub ls{

	my ($repertoire, $filtre) = @_;
	opendir(DOSSIER, $repertoire);

	if(!$filtre){

		$filtre = ".*";

	};	

	my @content = grep { !/^\.\.?$/ && -f "$repertoire/$_" && /$filtre/} readdir(DOSSIER);
	closedir(DOSSIER);

	return sort (map ({$repertoire.'/'.$_}  @content));
}


sub prefix_fasta_sequence{

	

	my ($file, $prefix) = @_;
	system("sed -i -e 's/>/>$prefix /;s/ /_/g' $file");
	#print "sed-i -e 's/>/> $prefix/' $file\n";

	return;
}

=head1 DEPENDENCIES

cap3

abyss-pe

=head1 INCOMPATIBILITIES

Fully compatible with any perl version

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 WRITTEN BY

=over 1

=item Sandy CONTRERAS (INRA), sandy.contreras-at-supagro.inra.fr

=back

=head1 EDITED BY

=over 1

=item Felix HOMA (INRA) felix.homa@cirad.fr

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
