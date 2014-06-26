#!/usr/bin/perl

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

annotator.pl -  make a blast of nucleic / proteic sequences on a custom nucleic / proteic databases list=head1 OPTIONS

=head1 DESCRIPTION

This script perform a blast of input sequences on the databases list using blast_custer.pl.
The no hit sequences of a blast are used as input for the following blast. 
If no database specified, blast will be performed on SwissProt, Trembl, NR and NT

=head1 SYNOPSIS / USAGE

annotator.pl -ni nucleic_sequences.fasta -pi proteic_sequences.fasta -o output_directory
[-nd|-nucleic_database	nucleic_database] 
[-pd|-proteic_database	-proteic_database] 
[-n|-name	analysis_name] 
[-q|-queue	sge_queue] \
[-e|evalue	evalue_cutoff] 
[-dust] 
[-max_target_seq maximum_number_of_aligned_sequences_to_keep] 
[-max_thread maximum_number_of_simultaneously_running_tasks] 
[-hold_jid sge_jobs_to_wait] 
[-num_seq_by_batch number_of_sequences_by_fasta_file] 
[-task task_to_execute] 
[-parse_deflines] 
[-task] 
[-advanced_option] 
[-clean remove repository] 

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Readonly;
use Getopt::Long;
use Pod::Usage;
#use Fatal qw/:void open close/;
use File::Path;

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

help "Not enough arguments" unless (@ARGV);

=pod

=head1 OPTIONS

=head2 required parameters

=over 11

=item B<[-ni]> <FASTA file>: 

Fasta file containing nucleic sequences

=item B<[-pi]> <FASTA file>: 

Fasta file containing proteic sequences

=item B<[-o]> <output directory>: 

Directory to store blast output

=back

=head2 optional arguments

=over 11

=item B<[-nd]> <nucleic database>: 

Database containing nucleic sequences

=item B<[-pd]> <proteic database>: 

Database containing proteic sequences

=item B<[-n|-name]> <analysis name>

Analysis name

=item B<[-q|-queue]> <sge queue>

SGE queue to use for blast, default bioinfo

=item B<[-e|-evalue]> <evalue cutoff>

Expectation value (E) threshold for saving hits, default 10

=item B<[-dust]>

Filter query sequence with DUST (Format: ’yes’,’level window linker’, or ’no’ to disable)

=item B<[-max_target_seq]>

Maximum number of aligned sequences to keep (<integer>), default 10

=item B<[-max_thread]> <number of threads>

Maximum number of simultaneously running tasks (<integer>), default 24

=item B<[-hold_jid]> <job to wait>

Wait for this job_id before processing the job_array (<integer>)

=item B<[-num_seq_by_batch]> <number of sequences>

Number of sequences by blast file (<integer>), default 16

=item B<[-parse_deflines]> <Boolean>

Should the query and subject defline(s) be parsed? 0 = no (Default), 1 = yes

=item B<[-advanced_option]> <options>

Advanced option, if you want to add more specific option (<string>), e.g : -advanced_option ’ -seg yes ’

=item B<[-clean]> <Boolean>

Remove repository. Default true

=item B<[-h|-help|-?]>

Prints a brief help message and exits.

=item B<[-man]>

Prints detailled help and script documentation.

=back

=cut


my($output, $help, $name, $man);
my @databases;
my @input;
my %analysis = ("nn" => "blastn", "np" => "blastx", "pn" => "tblastn", "pp" => "blastp");
my %blast_cluster;

$blast_cluster{'max_target_seq'}=10;
$blast_cluster{'num_seq_by_batch'}=16;

#$blast_cluster{'q'} = 'arcad.q';

GetOptions( 

	\%blast_cluster,		# Par défaut, si aucune variable n'est assignée, l'option est stockée dans le hash %blast_cluster

	'ni|nucleic_input=s'	=> sub{register($_[1], 'n', \@input)},
	'pi|proteic_input=s'	=> sub{register($_[1], 'p', \@input)},
	'o|output=s'			=> \$output,
	'nd|nucleic_database=s'	=> sub{register($_[1], 'n', \@databases)},
	'pd|proteic_database=s'	=> sub{register($_[1], 'p', \@databases)},
	'q|queue=s'				,
	'n|name=s'				=> \$name,
	'e|evalue=f'			,
	'dust=s'				,
	'max_target_seq=i'		,
	'max_thread=i'			,
	'hold_jid=i'			,
	'num_seq_by_batch=i'	,
	'task=s'				,
	'parse_deflines'		,
	'advanced_option=s'		,
	'clean=i'				,
	'help|h|?'				=> \$help,
	'man'					=> \$man
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

pod2usage(0) if(!@input || !$output);

if (!@databases){

	push (@databases, ['/SATA/bank/biomaj/uniprot/current/blast/uniprot_sprot.fasta', 'p'], ['/SATA/bank/biomaj/uniprot/current/blast/uniprot_trembl.fasta', 'p'], ['/SATA/bank/biomaj/nr/current/flat/nr', 'p'], ['/SATA/bank/biomaj/nt/current/flat/nt', 'n']);

}

my $step=1;

my $result_path = $output.'/'.$name.'/RESULT/';
mkpath($result_path);

my $blast_cluster_parameters="";

# Lecture des paramètres spécifiques à blast_cluster.pl
foreach my $param (keys(%blast_cluster)){ 

	$blast_cluster_parameters.=' -'.$param.' '.$blast_cluster{$param};

}


foreach my $database (@databases){ 

	# Pour chaque banque on créé un répertoire spécifique pour stocker les fichiers de blast et les résultats
	my ($database_path, $database_name, $database_extension)=($database->[0]=~/^([\.\/]+$|.+\/)*(\.[^\/\.]*|[^\/\.]*)(\..+)*$/);
	my $blast_path = $output.'/'.$name.'/BLAST_FILES/step'.$step.'_'.$database_name.'/';
	mkpath($blast_path);

	my @seq_to_blast;

	foreach my $file (@input){


		my ($file_path, $file_name, $file_extension)=($file->[0]=~/^([\.\/]+$|.+\/)*(\.[^\/\.]*|[^\/\.]*)(\..+)*$/);
		if(!defined($file_path)){$file_path="";print "$file_path\n";}
		if(!defined($file_extension)){$file_extension="";print "extension\n";}
		
		# Pour chaque fichier on réalise un blast sur la banque
		print 'blast_cluster.pl -input '.$file_path.$file_name.$file_extension.' -program '.$analysis{$file->[1].$database->[1]}.' -database '.$database->[0].' -output '.$file_name.'.xml -directory '.$blast_path.' -format 5 '.$blast_cluster_parameters."\n";

		system('blast_cluster.pl -input '.$file_path.$file_name.$file_extension.' -program '.$analysis{$file->[1].$database->[1]}.' -database '.$database->[0].' -output '.$file_name.'.xml -directory '.$blast_path.' -format 5 '.$blast_cluster_parameters); 

		
		my @xml_name = ($file_name=~/^(?:(.*?)(?:_no_hit))|(.*)/);

		# Pour chaque fichier on créé un xml qui contient l'ensemble des résultats de blast (sans les séquences no hit found)
		# et on récupère les séquences qui n'ont pas matché dans %no_hit_sequence
		my %no_hit_sequence=parseXML($blast_path.$file_name.'.xml', $result_path.join("", @xml_name).'_annotation.xml');

		# On génère un fichier fasta qui contient les séquences non matchées
		# et on ajoute ce fichier fasta à la liste des fichiers à blaster
		push(@seq_to_blast, [retrieve_No_Hit_Fasta_Sequence($file_path.$file_name.$file_extension, $blast_path.$file_name.'_no_hit_on_'.$database_name.'.fasta', \%no_hit_sequence), $file->[1]]);

	} 

	# La liste des fichiers à blaster (qui contiennent les séquences no hits) sera utilisé comme entrée pour le blast suivant
	@input = @seq_to_blast;
	$step++;
}


sub parseXML{

# Parse input xml to generate an output xml file containing only sequences with hit
# Return the list of sequences with no hit founds in a hash table

	my ($input_xml, $output_xml)=@_;
	
	my %no_hit_ids;

	open (INPUT_XML, '< '.$input_xml) or die ("\nCannot read xml file: $!\n");
	open (OUTPUT_XML, '>> '.$output_xml) or die ("\nCannot create xml file: $!\n");

	while(my $line=<INPUT_XML>){

		my $hit="true";
		my $line_to_print="";
		my $id_seq="";

		if ($line =~ /<Iteration>/){


			do{
	
				$line_to_print.= $line;

				if($line =~ /No hits found/){$hit="false";}

				if($line =~ /<Iteration_query-def>\s*(.*)\s*<\/Iteration_query-def>/){$id_seq=$1}
					
				$line=	<INPUT_XML>;

			}while ($line!~/<\/Iteration>/);

			# Si la séquence a matché, on l'inscrit dans le fichier XML de résultat
			if($hit eq "true"){print OUTPUT_XML $line_to_print.$line}

			# Sinon, on ajoute l'identifiant de la séquence dans une table de hash
			else{$no_hit_ids{$id_seq}=1;}

		}

		

	}

	close(INPUT_XML);
	close(OUTPUT_XML);

	return %no_hit_ids;

}

sub retrieve_No_Hit_Fasta_Sequence{

# Parse a FASTA file and use and use an hash table of id to to generate an output FASTA file containing only no hit sequences

	my ($input_fasta, $output_fasta, $id_list)=@_;
	open (INPUT_FASTA, '< '.$input_fasta) or die ("\nCannot read fasta file: $!\n");
	open (OUTPUT_FASTA, '>> '.$output_fasta) or die ("\nCannot create fasta file: $!\n");

	while (my $line=<INPUT_FASTA>){

		my $id_sequence="no";
		my $line_to_print="";

		if ($line=~/>\s*(.*)\s*/){$id_sequence=$1;};

		while($line && $id_list->{$id_sequence}){

			$line_to_print.=$line;
			$line=<INPUT_FASTA>;
			if ($line && $line=~/>\s*(.*)\s*/){$id_sequence=$1}

		}

		print OUTPUT_FASTA $line_to_print;

	}

	close(INPUT_FASTA);
	close(OUTPUT_FASTA);

	return $output_fasta;

}

sub register{

	my ($parameter, $type, $analysis)=@_;
	push (@$analysis, [$parameter, $type]);

}

=head1 DEPENDENCIES

blast_cluster.pl

=head1 INCOMPATIBILITIES

Fully compatble with any perl version

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 2

=item Sandy CONTRERAS (INRA), sandy.contreras-at-supagro.inra.fr

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
