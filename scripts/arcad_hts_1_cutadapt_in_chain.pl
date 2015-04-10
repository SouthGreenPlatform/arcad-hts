#!/bin/env perl

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

=head1 NAME

cutadapt_in_chain.pl - Launch Cutadapt on fastq files

=head1 DESCRIPTION

Launch Cutadapt on all fastq of a directory and its arborescence if requested. A global or a specific adapter set of sequences to trim can be provided by the user.

=head1 SYNOPSIS / USAGE

cutadapt_in_chain.pl	-i input_folder \
[-o output_folder] \
[-zip|nozip] \
[-sub|-nosub] \
[-reverse|-noreverse] \
[-a global_file_for_a_adapter] \
[-g global_file_for_a_adapter] \
[-b global_file_for_b_adapter] \
[-sa specific_file_for_a_adapter] \
[-sg specific_file_for_b_adapter] \
[-sb specific_file_for_b_adapter] \
[-q quality_level] \
[-v overlap] \
[-e error_rate] \
[-u number_of_bases_to_hard_trim] \
[-queue SGE_queue]

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Readonly;
use Getopt::Long;
use Pod::Usage;
use Fatal qw/:void open close/;
use File::Basename;

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

=over 1

=item B<[-i]> ([input_folder]): 

The directory containing the fastq files in it or in its arborescence

=back

=head2 optional arguments

=over 11

=item B<[-o]> ([output_folder]) [.]

If provided, Cutadapt outputs will be gather in this directory with file having the same name as the input. 
Default behavior is to have the Cutadapt outputs in the same directory as the fastq file.


=item B<[-sub|nosub]> ([flag to check subdirectory]) [-sub]

-sub will allow to launch Cutadapt on fastq files present in subdirectories
-nosub disable this possibility
Default: -sub is on

=item B<[-zip|nozip]> ([flag to check subdirectory]) [-zip]

-zip will compress output files with gzip
-nozip will output fastq files
Default: -zip is on

=item B<[-reverse|-noreverse]> ([flag to reverse adapter sequences]) [-reverse]

-reverse will make cutadapt searching also for reverse complement of adapter sequences
-noreverse disable this possibility
Default -reverse is on


=item B<[-a]> adapter file []

This option defines the adapters to search at the end of the read in every fastq file.
Adpater file must be in fasta format.
At least one of -a/-b/-sa/-sb option must be specified

=item B<[-g]> adapter file

This option defines the adapters to search at the beginning of the read in every fastq file.
Adpater file must be in fasta format.
At least one of -a/-b/-sa/-sb option must be specified

=item B<[-b]> adapter file

This option defines the adapters to search everywhere in the read in every fastq file.
Adpater file must be in fasta format.
At least one of -a/-b/-sa/-sb option must be specified


=item B<[-sa]> adapter file

This option defines the adapters to search at the end of the read. Each directory containing fastq file must contain an adapter file with this name.
Adpater file must be in fasta format.
At least one of -a/-b/-sa/-sb option must be specified

=item B<[-sg]> adapter file

This option defines the adapters to search at the beginning of the read. Each directory containing fastq file must contain an adapter file with this name.
Adpater file must be in fasta format.
At least one of -a/-b/-sa/-sb option must be specified

=item B<[-sb]> adapter file

This option defines the adapters to search everywhere in the read. Each directory containing fastq file must contain an adapter file with this name.
Adpater file must be in fasta format.
At least one of -a/-b/-sa/-sb option must be specified


=item B<[-q]> quality threshold [20]

The quality threshold to pass to cutadapt for trimming the end of the reads.
Default: 20


=item B<[-v]> minimum overlap [7]

Minimum overlap to detect the adapter
Default: 7


=item B<[-e]> error ratio [0.1]

Error ratio allowed to detect the adapter.
Default: 0.1


=item B<[-m]> minimum read size [20]

Minimum read size to keep. Reads shorter than this value will be discarded
Default: 20

=item B<[-u]> bases to remove [0]

The number of bases to remove from the read. A positive number will remove bases from the
beginning of the read and a negative number from the end of the read. Can be specified 
once for each sign 
Default: No trimming

=item B<[-queue]> SGE queue [bioinfo.q]

The SGE queue to launch cutadapt
Default: bioinfo.q

=back

=cut

my $version="";
my $CUTADAPT_PATH = &$Softwares::CUTADAPT_PATH or confess("$!");
my $MAX_PARALLEL=12;

#options processing 
my ($man, $help, $debug, $input, $output, $subdirectory, $reverse, $a, $g, $b, $sb, $sg, $sa, $quality, $overlap, $minsize, $error, $queue, @cut, $zipoutput);
$queue = "normal.q";
$subdirectory=1;
$zipoutput=1;
$reverse=1;
$quality=20;
$overlap=7;
$error=0.1;
$minsize=20;
my $ra_files;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
#	"debug:i"     => \$debug,
	"input|i=s"   => \$input,
	"output|o=s"  => \$output,
	"zip!"        => \$zipoutput,
	"sub!"        => \$subdirectory,
	"reverse!"    => \$reverse,
	"a=s"         => \$a,
	"g=s"         => \$g,
	"b=s"         => \$b,
	"sb=s"        => \$sb,
	"sg=s"        => \$sg,
	"sa=s"        => \$sa,
	"quality|q=i" => \$quality,
	"overlap|v=i" => \$overlap,
	"minsize|m=i" => \$minsize,
	"error|e=f"   => \$error,
	"cut|u=s"        => \@cut,	
	"queue=s"     => \$queue
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

pod2usage(0) if(!(defined($a)||defined($b)||defined($sa)||defined($sb)||defined($sg)||defined($g)));

if(defined($output) && -e $output && !(-d $output))
{
	warn "$output already exists, but it's not a directory";
	exit(0);
}
mkdir $output if(defined($output) && !(-e $output));

if( !-e $input ){
	help "$! $input";
}
elsif( -f $input ){
	push @$ra_files, $input;
}
elsif( -d $input )
{
	$ra_files = Files->getFiles( $input, '*.fastq', $subdirectory );
	$ra_files = Files->getFiles( $input, '*.fq', $subdirectory );
	$ra_files = Files->getFiles( $input, '*.fastq.gz', $subdirectory );
	$ra_files = Files->getFiles( $input, '*.fq.gz', $subdirectory );
	help "No FASTQ file(s) in input" if( 0 >= scalar @$ra_files );
}
else{
	help "$! $input";
}
if(scalar(@cut) >2)
{
	help "You provided too many times the cut/u options";
}
elsif(scalar(@cut) == 2)
{
	if(($cut[0] > 0 && $cut[1] > 0) || ($cut[0] < 0 && $cut[1] < 0))
	{
		help "You provided too many times the cut/u options with the same sign";
	}
	
}
my $adapters="";
if(defined($a))
{
	my $adaptersTmp = &readAdapters( "$a" );
	$adaptersTmp =~ s/\// -a /g;
	$adapters .= ' -a ' . $adaptersTmp . ' ';
}

if(defined($b))
{
	my $adaptersTmp = &readAdapters( "$b" );
	$adaptersTmp =~ s/\// -b /g;
	$adapters .= ' -b ' . $adaptersTmp;
}

if(defined($g))
{
	my $adaptersTmp = &readAdapters( "$g" );
	$adaptersTmp =~ s/\// -g /g;
	$adapters .= ' -g ' . $adaptersTmp;
}


my $i=0;
foreach my $fastq (@$ra_files)
{
	
	my($file, $path, $ext) = fileparse($fastq, qr/\.[^.]*[\.gz]?/);
	$file .= $ext;

	my $spe_adapters="";
	if(defined($sa))
	{
		my $adaptersTmp = &readAdapters( "$path/$sa" );
		$adaptersTmp =~ s/\// -a /g;
		$spe_adapters .= ' -a '.$adaptersTmp;
	}
	if(defined($sb))
	{

	        my $adaptersTmp = &readAdapters( "$path/$sb" );
		$adaptersTmp =~ s/\// -b /g;
		$spe_adapters .= ' -b '.$adaptersTmp;
	}
	if(defined($sg))
	{
		my $adaptersTmp = &readAdapters( "$path/$sg" );
		$adaptersTmp =~ s/\// -g /g;
		$spe_adapters .= ' -g '.$adaptersTmp;
	}
	my ($output_path, $out_file);
	
	if($zipoutput && $ext !~ /\.gz$/)
	{
		$out_file = $file.".gz";
	}
	elsif(! $zipoutput && $ext =~ /\.gz$/)
	{
		$out_file = substr($file, -3);
	}
	else
	{
		$out_file = $file;
	}
	
	if(defined($output)){
		$output_path="$output/$out_file";
	}
	else
	{
		$out_file=~s/\.(f[ast]?q)/\.cutadapt\.$1/;
		$output_path="$path/$out_file";
	}
	
	my $cut="";
	foreach my $totrim (@cut)
	{
		$cut.= " -u $totrim";
	}
	
	print "OK\n";
	if($i>=$MAX_PARALLEL)
	{
		my $ok=0;
		while(!$ok)
		{
			if((my $cnt =`qstat|grep -c "cut_arcad"`) < $MAX_PARALLEL){
				$ok=1;
			}
			else{
				sleep(5);
			}
		}
	}
	
	print STDERR ("qsub -q $queue -b yes -V -cwd -N cut_arcad $CUTADAPT_PATH $adapters $spe_adapters $cut -o $output_path -q $quality -O $overlap -m $minsize -e $error $path/$file");
	system("qsub -q $queue -b yes -V -cwd -N cut_arcad $CUTADAPT_PATH $adapters $spe_adapters $cut -o $output_path -q $quality -O $overlap -m $minsize -e $error $path/$file");
	$i++;
}

### This function read adapters file and concat all adapters.
### By default seperation caracter is '/' it can be changed.
sub readAdapters
{
	my $filename = shift;
	my $sep = '/';
	$sep = shift if( @_ );
	
	my $adapt = '';
	open(my $file_handle, $filename) or confess("$!");
		while(my $seq=<$file_handle>)
		{
			chomp $seq;
			$seq =~ s/\r//;
			next if $seq=~m/^>/;#Format fasta
			$adapt.="$seq$sep";
			if($reverse)
			{
				my $revcomp = reverse($seq);
				$revcomp =~ tr/ACGTacgt/TGCAtgca/;
				$adapt.="$revcomp$sep";
			}
		}
		$adapt =~ s/$sep$//;
	close $file_handle;
	return $adapt;
}

=head1 DEPENDENCIES

Cutadapt

=head1 INCOMPATIBILITIES

Fully compatble with any perl version

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 2

=item Gautier SARAH (INRA), gautier.sarah-at-supagro.inra.fr

=item Francois SABOT (IRD), francois.sabot-at-ird.fr

=item Felix Homa (INRA), felix.homa@cirad.fr

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

=head1 TODO

Print warning message if input directory is same than output directory.
Because input data are overwritten.

=cut
