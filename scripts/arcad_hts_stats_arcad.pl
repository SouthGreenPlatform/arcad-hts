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

Stats_arcad - Compute statistics from the different arcad analysis

=head1 DESCRIPTION

This script makes a statistical summary over a NGS data analysis in a file

=head1 SYNOPSIS / USAGE

Stats_arcad.pl -b bam_file -c coverage_file -p prot4est -f fasta_file -blast blast_file (as much as needed) -o output_file obsolete
Stats_arcad.pl -b bam_file -p prot4est -f fasta_file -blast blast_file (as much as needed) -o output_file 

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Readonly;
use Getopt::Long;
use Pod::Usage;
use Fatal qw/:void open close/;
use Cwd qw(cwd abs_path);
use Tie::File;

use lib '/NAS/arcad_data/Softs/tags';
use Modules::Config::Softwares;

=head1 OPTIONS

=head2 required parameters

=over 5

=item B<[-b]> ([a BAM file]): 

The BAM alignment file

=item B<[-f]> ([a fasta file containing the reference sequences]):

=item B<[-o]> ([output file])
File were the statistics will be written

=head2 optional arguments

=over

=item B<[-p]> ([the prot_main.psql file from prot4est]):

=item B<[-blast]> ([a blast output file at XML format])
This option can be set as many times as needed

=back

=cut

my $version="";

unless (@ARGV)
{
  pod2usage(0);
  exit();
}

my $SAMTOOLS_PATH = &$Softwares::SAMTOOLS_PATH or confess("$!"); # "samtools";

#options processing
my ($man, $help, $debug, $bam, $coverage, $prot4est, $reference, @blast, $output);

GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "bam|b=s"    => \$bam,
           "p4e|p=s"    => \$prot4est,
           "fasta|f=s"  => \$reference,
           "blast=s"    => \@blast,
           "output|o=s" => \$output 
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

#create the tmp directory in the same directory as the output
my $tmp = $output; 
$tmp =~ s/[^\/]+$/tmp/;
my $log="$tmp/../log_stats";
mkdir($tmp);

$tmp = abs_path($tmp);
$log = abs_path($log);

system("echo -e \"Start\n\" >$log && date >>$log" );
#Checks the presence of the bam index file. If it does not exists, it is created.
my $bam_index=$bam;
$bam_index =~ s/bam$/bai/;
if( ! -e $bam_index )
{
  system("$SAMTOOLS_PATH index $bam $bam_index"); 
}

$bam = abs_path($bam);

system("echo -e \"extracting sequence length\n\" >>$log && date >>$log" );
extract_sequence_length($bam);



system("echo -e \"extracting mapping length\n\" >>$log && date >>$log" );
extract_mapping_length($bam);

if(defined($prot4est) && -e $prot4est)
{
  system("echo -e \"extracting translation informations\n\" >>$log && date >>$log" );
  extract_p4e_information($prot4est);
}

extract_individual_FPKM($bam, $tmp, $log);

concatenate_statistics();

=head1 FUNCTIONS

=head2 extract_sequence_length

=head3 args

=item bam file - Bam file to be processed

=head3 return

=item void

=cut

sub extract_sequence_length
{
  my $bam = $_[0];
  system("$SAMTOOLS_PATH idxstats $bam | awk '{print \$1\"\t\"\$2}' >$tmp/length.stats");
}

=head2 extract_mapping_length

=head3 args

=item bam file - Bam file to be processed

=head3 return

=item void

=cut

sub extract_mapping_length
{
  my $bam = $_[0];
  system("$SAMTOOLS_PATH depth $bam | awk '{if( \$3 > 0 ){ counter[\$1]++ }} END{for(i in counter){print i\"\t\"counter[i] } }' > $tmp/mapping_length.stats");
}

=head2 extract_p4e_information

=head3 args

=item prot4EST output file - prot4EST output file to be processed

=head3 return

=item void

=cut

sub extract_p4e_information
{
  my $prot4est=$_[0];
  system("awk -F \",\" '{start=0;if(\$5>=0){start=\$5}else{start=\$6};end=0;if(\$9>=0){end=\$9}else{end=\$8};if(\$7>0){print \$2\"\t\"\$3\"\t\"start\"\t\"end\"\t\"\$7}else{print \$2\"\t\"\$3\"\t\"end\"\t\"start\"\t\"\$7}}' $prot4est > $tmp/prot.stats");
}

=head2 extract_individual_FPKM

=head3 args

=item bam file - Bam file to be processed

=item temporary folder - Folder in which global bam will be be decomposed in many bams

=item log file - File in which all error and residual output will be log

=head3 return

=item void

=cut

sub extract_individual_FPKM
{
  my ($global_bam, $tmp, $log) = @_;
  
  my $base_dir = cwd;
  #Extraction des nom des readgroup
  system("$SAMTOOLS_PATH view -H $global_bam |grep -v bwa |awk '/\@RG/{sub(\"ID:\",\"\",\$2);print \$2}' >$tmp/RG");
  #Creation d'un fichier BAM par readgroup
  my $bam_dir = "$tmp/bam";
  mkdir("$bam_dir");
  chdir "$bam_dir";
  
  my $file_handle;
  if ( open($file_handle, "$tmp/RG") )
  {
    while(<$file_handle>)
    {
      system("echo -e \"Calculating RPKM for $_\n\" >>$log && date >>$log" );
      chomp();
      system("$SAMTOOLS_PATH view -b -r $_ -o $bam_dir/$_.bam $global_bam");
      system("$SAMTOOLS_PATH index $bam_dir/$_.bam $bam_dir/$_.bai");
      system("$SAMTOOLS_PATH idxstats $_.bam >$_.stat");
      system("awk '{total+=\$3;mapped[\$1]=\$3;len[\$1]=\$2} END{for(x in mapped){if(len[x]!=0){rpkm=(mapped[x])/(len[x]/1000*total/1000000);print x\"\t\"rpkm}}}' $bam_dir/$_.stat > $tmp/$_.rpkm.stats");
    }
  }
  close($file_handle);
  chdir "$base_dir";
}

=head2 concatenate_statistics

=head3 args

=item void

=head3 return

=item void

=cut

sub concatenate_statistics
{
	my $header="Contig\tlength\tmapping length\tp4e method\tCDS start\tCDS end\tframe";
	# Ouverture du fichier de sorti
	my %h_all_info = ();
	my @a_all_rpkm_header = ('LEN', 'COV', 'METH', 'START', 'END', 'FRAME');
	
	if( -e "$tmp/length.stats" )
	{
		open( my $TXT, "$tmp/length.stats") or confess("$!");
			while(my $line = <$TXT>)
			{
				chomp($line);
				if( $line =~ /(\S+)\s(\d+)/ )
				{
					my($nom, $taille) = ($1, $2);
					${$h_all_info{$nom}{LEN}} = $taille;
				}
			}
		close($TXT);
	}
	
	if( -e "$tmp/mapping_length.stats" )
	{
		open( my $TXT, "$tmp/mapping_length.stats") or confess("$!");
			while(my $line = <$TXT>)
			{
				chomp($line);
				if( $line =~ /(\S+)\s(\d+)/ )
				{
					my($nom, $coverage) = ($1, $2);
					${$h_all_info{$nom}{COV}} = $coverage;
				}
			}
		close($TXT);
	}
	
	if(-e "$tmp/prot.stats")
	{
		open( my $TXT, "$tmp/prot.stats") or confess("$!");
			while(my $line = <$TXT>)
			{
				chomp($line);
				if( $line =~ /(\S+)\s(\w+)\s(\d+)\s(\d+)\s([\-\+\d+]+)/ )
				{
					my($nom, $method, $cds_start, $cds_end, $frame) = ($1, $2, $3, $4, $5);
					${$h_all_info{$nom}{METH}}  = $method;
					${$h_all_info{$nom}{START}} = $cds_start;
					${$h_all_info{$nom}{END}}   = $cds_end;
					${$h_all_info{$nom}{FRAME}} = $frame;
				}
			}
		close($TXT);
	}
	
	if (open(my $file_handle, "$tmp/RG") )
	{
		while(<$file_handle>)
		{
			chomp();
			my $file = $_;
			if(-e "$tmp/$file.rpkm.stats")
			{
				open( my $TXT, "$tmp/$file.rpkm.stats") or confess("$!");
					$header .= "\tRPKM $file";
					push @a_all_rpkm_header, "RPKM $file";
					while(my $line = <$TXT>)
					{
						chomp($line);
						if( $line =~ /(\S+)\s([\+\-\.\d+]+)/ )
						{
							my($nom, $rpkm) = ($1, $2);
							if( !exists $h_all_info{$nom} ){
								$h_all_info{$nom} = {"RPKM $file" => $rpkm};
							}
							else{
								${$h_all_info{$nom}{"RPKM $file"}} = $rpkm;
							}
						}
					}
				close $TXT;
			}
		}
	}
	
	open(my $TXT, ">$output") or confess("$!");
	print $TXT $header, "\n";
	foreach my $key (sort {$a cmp $b} keys %h_all_info)
	{
		print $TXT "$key\t";
		foreach my $subKey ( @a_all_rpkm_header )
		{
			if( exists $h_all_info{$key}{$subKey} ){
				print $TXT ${$h_all_info{$key}{$subKey}} . "\t";
			}
			else{
				print $TXT ".\t";
			}
		}
		print $TXT "\n";
	}
	close $TXT
}

=head1 DEPENDENCIES

Unix system only
Samtools

=head1 INCOMPATIBILITIES

Windows, Mac

=head1 BUGS AND LIMITATIONS

Void

=head1 AUTHORS

=over 2

=item Gautier SARAH

=item gautier.sarah-at-supagro.inra.fr

=item Felix Homa

=item felix.Homa-at-cirad.fr

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
