#! /bin/env perl
#  
#  Copyright 2014 INRA-CIRAD-IRD
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

=head1 Map Fastq Files


Mapping_Arcad - Map Fastq Files on a reference and produce one merged BAM file

=head1 SYNOPSIS

mapping_arcad.pl -i conf_file -o bam_file -r reference_fasta_file -keep_intermediate -rmdup

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Readonly;
use Getopt::Long;
use Pod::Usage;
use Fatal qw/:void open close/;
use File::Path qw (rmtree);
use File::Basename qw (fileparse);

use lib '/NAS/arcad_data/Softs/tags';
use Modules::Config::Softwares;

=pod

=head1 OPTIONS


=head2 Parameters

=over 4

=item B<[-i]> ([a conf file]): 

This file must be tabulated. 
First column must contain a path to a fastq file
Second column must contain the RG ID associated to the fastq file
Third column must contain "single" or "paired"
Fourth column if data are paired must contain the associated fastq file of the first one

=item B<[-o]> ([a BAM File]): 

The output BAM file of the whole process

=item B<[-r]> ([a fasta file containing the reference sequences])

A Fasta file containing the reference.

=item B<[-q]> ([a queue])

A valid SGE queue.
Default is arcad.q

=item B<[-keep_intermediate/nokeep_intermediate]> ([a flag ])

Whether or not keep intermediate files
By default keep_intermediate is OFF

=item B<[-rmdup/normdup]> ([a flag])

Whether or not process rmdup on the alignments (As this must be done before merging them)
By default rmdup is ON

=item B<[-mapper]> ([mapper name])

Can be one of these : bwa, bwa_mem, gem

=cut

my $version="";

pod2usage(0) unless (@ARGV);

#options processing
my ($man, $help, $debug, $conf_file, $output, $reference, $keep_intermediate, $rmdup, $queue, $mapper );
$keep_intermediate = 0;
$rmdup = 1;
$queue="arcad.q";
$mapper='bwa';

GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "keep_intermediate|k!" => \$keep_intermediate,
           "rmdup!"     => \$rmdup,
           "input|conf|i=s" => \$conf_file, 
           "reference|r=s"     => \$reference,
           "queue|q=s"   => \$queue,
           "output|o=s"  => \$output,
	   "mapper:s"       => \$mapper
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
pod2usage({	-msg     => "Please set valid value for options -r and -i"
		-exitval => 0
	}
) if(!(defined($reference) && defined($conf_file) && defined($output)));
	
pod2usage({	-msg     => "Cannot find the reference file reference and/or config file and/or you did not set the output\nPlease set valid value for options -r and -i",
		-exitval => 0
	}
)if( !(-e $reference && -e $conf_file) );

my %MAPPERS = (bwa => \&map_bwa, bwa_mem => \&map_bwa_mem );
my $PICARD_TOOLS_DIRECTORY=&$Softwares::PICARD_TOOLS_DIRECTORY;
my $JAVA_PATH=&$Softwares::JAVA_PATH;

my $HEADER = "#! $ENV{SHELL}\n#\$ -q $queue\n";
$HEADER .= <<'HEADER';
#$ -cwd
#$ -V
#$ -b yes
#$ -pe parallel_smp 4
#$ -N Map_Arc

HEADER

#Parse conf file and process the different mapping
my ($fi, $tmp_directory, $out_extention) = fileparse($output, qr/\.[^.]*/);
$tmp_directory .= "${fi}_tmp";
mkdir($tmp_directory);
#my $conf_handle;
my %mapping_files;
if (open(my $conf_handle, $conf_file))
{
	while(<$conf_handle>)
	{
		chomp();
		my @splitted = split(/\t/, $_);
		if(! -e $splitted[0])
		{
			print "Cannot find the fastq file ".$splitted[0];
			exit();
		}
		if(scalar(@splitted == 4) && ! -e $splitted[3])
		{
			print "Cannot find the fastq file ".$splitted[3];
			exit();
		}
		my $bam_output = "$splitted[1].$splitted[2].bam";
		my $i=1;
		while(exists($mapping_files{"$tmp_directory/$bam_output"}))
		{
			$bam_output= "$splitted[1].$splitted[2]_$i.bam";
			$i++;
		}
		$mapping_files{"$tmp_directory/$bam_output"}++;
		
		if(scalar(@splitted==3))
		{
			$MAPPERS{$mapper}->($splitted[0], $reference, $bam_output, $splitted[1], $rmdup,$tmp_directory);
		}
		elsif(scalar(@splitted==4))
		{
			$MAPPERS{$mapper}->($splitted[0], $reference, $bam_output, $splitted[1], $rmdup, $tmp_directory, $splitted[3]);
		}
	}
	close($conf_handle);
	
	#Wait for all mapping to end
	sleep(10) while( (my $cnt =`qstat|grep -c "Map_Arc"`) > 0 );
	sleep(30);
	
	#merge the different bam_files
	merge_mappings(\%mapping_files, $output) if( $out_extention eq '.bam' || $out_extention eq '.sam' );
	rmtree($tmp_directory) if(!$keep_intermediate);
}
else
{
	print "Cannot open $conf_file";
	exit();
}

exit;

####################
#
# SUB PROGRAMS
#
####################

#Create a sh file and launch the single end mapping

sub map_bwa
{
	my ($forward, $reference, $bam_output,$rg,$rmdup,$directory, $reverse)=@_;
	my($pattern, undef, undef) = fileparse($bam_output, qr/\.[^.]*/);
	my $command_file="$pattern.sh";
	
	my $BWA_PATH =&$Softwares::BWA_PATH;
	my $BWA_ALN=$BWA_PATH . " aln -n 3 -t 4";
	my $BWA_SAMSE=$BWA_PATH . " samse";
	my $BWA_SAMPE=$BWA_PATH . " sampe -a 1000 ";
	# Cette ligne bash arrete l'execution d'un programme qui a foire ‡ une etape
	my $TEST_BASH = 'if [[ $? -gt 0 ]] ; then echo "program killed"; exit $?; fi';
	
	if(! -e "$reference.ann")
	{
		system("qsub -q $queue -sync yes -b yes -N bwa_index bwa index -a is $reference");
		#print "[BWA] Indexing $reference OK\n" if( $debug );
	}

	if(open(my $com_handle,">$directory/$command_file"))
	{
		print $com_handle $HEADER;
		print $com_handle $BWA_ALN." $reference ".$forward." > $directory/$pattern.forward.sai\n\n";
		print $com_handle "$TEST_BASH\n\n";
		#print $com_handle "echo \"[BWA] Indexing $forward OK\n\"" if( $debug );
		if( defined $reverse && -e $reverse )
		{
			print $com_handle $BWA_ALN." $reference ".$reverse." > $directory/$pattern.reverse.sai\n\n";
			print $com_handle "$TEST_BASH\n\n";
			#print $com_handle "echo \"[BWA] Indexing $reverse OK\n\"" if( $debug );
			print $com_handle $BWA_SAMPE." -r \'\@RG\tID:".$rg."\tPL:ILLUMINA\tSM:".$rg."\' $reference $directory/$pattern.forward.sai $directory/$pattern.reverse.sai $forward $reverse > $directory/$pattern.sam\n\n";
			print $com_handle "$TEST_BASH\n\n";
			#print $com_handle "echo \"[BWA] Paired end mapping OK\n\n\"" if( $debug );
		}else{
			print $com_handle $BWA_SAMSE." -r \'\@RG\tID:".$rg."\tPL:ILLUMINA\tSM:".$rg."\' $reference $directory/$pattern.forward.sai ".$forward." > $directory/$pattern.sam\n\n";
			print $com_handle "$TEST_BASH\n\n";
			#print $com_handle "echo \"[BWA] Single end mapping OK\n\n\"" if( $debug );
		}
		print $com_handle $JAVA_PATH." -jar $PICARD_TOOLS_DIRECTORY/SortSam.jar I=$directory/$pattern.sam O=$directory/$bam_output SO=coordinate VALIDATION_STRINGENCY=SILENT\n\n";
		print $com_handle "$TEST_BASH\n\n";
		#print $com_handle "echo \"[PICARD] Sort $directory/$pattern.sam by coordinate OK\n\n\"" if( $debug );
		if($rmdup)
		{
			print $com_handle "$JAVA_PATH -Xmx2g -jar $PICARD_TOOLS_DIRECTORY/MarkDuplicates.jar I=$directory/$bam_output O=$directory/$bam_output.tmp VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE M=$bam_output.metrics \n";
			print $com_handle "$TEST_BASH\n\n";
			print $com_handle "mv -f $directory/$bam_output.tmp $directory/$bam_output\n";
			print $com_handle "$TEST_BASH\n\n";
			#print $com_handle "echo \"[PICARD] Rmdup $OK\n\n\"" if( $debug );
		}
		close($com_handle);
		chmod(0755 , "$directory/$command_file");
		system("qsub  $directory/$command_file");
	}
}

sub map_bwa_mem
{
	my ($forward, $reference, $bam_output, $rg, $rmdup, $directory, $reverse)=@_;
	my($pattern, undef, undef) = fileparse($bam_output, qr/\.[^.]*/);
	my $command_file="$pattern.sh";
	my $BWA_PATH =&$Softwares::BWA_PATH;
	# Cette ligne bash arrete l'execution d'un programme qui a foire ‡ une etape
	my $TEST_BASH = 'if [[ $? -gt 0 ]] ; then echo "program killed"; exit $?; fi';
	
	if(! -e "$reference.ann")
	{
		system("qsub -q $queue -sync yes -b yes -N bwa_index bwa index -a is $reference");
		#print "[BWA] Indexing $reference OK\n" if( $debug );
	}

	if(open(my $com_handle,">$directory/$command_file"))
	{
		print $com_handle "$HEADER\n";
		if( defined $reverse && -e $reverse )
		{
			print $com_handle $BWA_PATH . " mem -C -M -R '\@RG\tID:$rg\tPL:ILLUMINA\tSM:$rg' $reference $forward $reverse > $directory/$pattern.sam\n\n";
			print $com_handle "$TEST_BASH\n\n";
			#print $com_handle "echo \"[BWA] Paired end mapping OK\n\n\"" if( $debug );
		}else{
			print $com_handle $BWA_PATH . " mem -C -M -R '\@RG\tID:$rg\tPL:ILLUMINA\tSM:$rg' $reference $forward > $directory/$pattern.sam\n\n";
			print $com_handle "$TEST_BASH\n\n";
			#print $com_handle "echo \"[BWA] Single end mapping OK\n\n\"" if( $debug );
		}
		print $com_handle $JAVA_PATH." -jar $PICARD_TOOLS_DIRECTORY/SortSam.jar I=$directory/$pattern.sam O=$directory/$bam_output SO=coordinate VALIDATION_STRINGENCY=SILENT\n\n";
		print $com_handle "$TEST_BASH\n\n";
		#print $com_handle "echo \"[PICARD] Sort $directory/$pattern.sam by coordinate OK\n\n\"" if( $debug );
		if($rmdup)
		{
			print $com_handle "$JAVA_PATH -Xmx2g -jar $PICARD_TOOLS_DIRECTORY/MarkDuplicates.jar I=$directory/$bam_output O=$directory/$bam_output.tmp VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE M=$bam_output.metrics \n";
			print $com_handle "$TEST_BASH\n\n";
			print $com_handle "mv -f $directory/$bam_output.tmp $directory/$bam_output\n";
			print $com_handle "$TEST_BASH\n\n";
			#print $com_handle "echo \"[PICARD] Rmdup $OK\n\n\"" if( $debug );
		}
		close($com_handle);
		chmod(0755 , "$directory/$command_file");
		system("qsub  $directory/$command_file");
	}
}


#Merge all the generated BAM files
sub merge_mappings
{
	my ($mappings,$output) = @_;
	my $input="";
	$input.="I=$_ " for( keys(%$mappings) );
	my $command = $JAVA_PATH." -jar $PICARD_TOOLS_DIRECTORY/MergeSamFiles.jar $input O=$output MSD=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SO=coordinate";
	system("qsub -q $queue -b yes -V -cwd -sync yes -N Merge $command");
}
	
=pod

=head1 AUTHORS

Gautier SARAH (INRA), gautier.sarah@supagro.inra.fr
Fran√ßois SABOT (IRD), francois.sabot@ird.fr

=head1 VERSION
Version 1.0 

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
