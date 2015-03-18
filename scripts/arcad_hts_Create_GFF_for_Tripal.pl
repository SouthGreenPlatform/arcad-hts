#! /bin/env perl

#  
#  Copyright 2014 INRA
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

=head1 Create GFF for Tripal

Create a gff3 file for tripal 

=head1 SYNOPSIS

Create_GFF_for_Tripal.pl -f fasta_file -b blast2go  -xn blast_name -xf blast_xml_file -s stat_file -p prefix

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Readonly;
use Getopt::Long;
use Pod::Usage;
use Fatal qw/:void open close/;
use File::Path qw (rmtree);

use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::SearchIO; 
use Bio::Tools::GFF;

use Data::Dumper;
=pod

=head1 OPTIONS


=head2 Parameters

=over 5

=item B<[-p]> ([a string]): 

A prefix that will be added to all sequences name

=item B<[-f]> ([a fasta file]): 

The fasta file containing contigs

=item B<[-b]> ([a blast2GO annot file]): 

The annot file from blast2GO

=item B<[-xn]> ([a blast name]): 

A blast name. -xn and -xf must be provided in the same order

=item B<[-xf]> ([a blast XML file]): 

The blast xml file. -xn and -xf must be provided in the same order

=item B<[-s]> ([an arcad stat file]): 

The tabulated arcad stat file.

=item B<[-o]> ([gff ouptut]): 

The gff output file

=cut

my $version="";

unless (@ARGV)
{
  pod2usage(0);
}

my ($man, $help, $debug, $fasta, $prefix, $annot, $stat, @blast_name, @blast_file, $output);

$prefix="";

# parse options and print usage if there is a syntax error.
#--- see doc at http://perldoc.perl.org/Getopt/Long.html
GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "prefix|p=s" =>\$prefix,
           "fasta|f=s" =>\$fasta,
           "b2g|b=s" => \$annot, 
           "stat|s=s" => \$stat,
           "blast_name|xn=s" => \@blast_name, 
           "blast_file|xf=s"  => \@blast_file,
           "output|o=s"  => \$output,
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

if(! -e $fasta) 
{
  print("Cannot find $fasta");
  pod2usage(0);
}

if(defined($annot) && ! -e $annot)
{
  print("Cannot find $annot");
  pod2usage(0);
}

if(defined($stat) && ! -e $stat)
{
  print("Cannot find $stat");
  pod2usage(0);
}

if(!defined($output))
{
  print("You did not precise an output file");
  pod2usage(0);
}

if(scalar(@blast_name)!=scalar(@blast_file))
{
  my $nb_name = scalar(@blast_name);
  my $nb_file = scalar(@blast_file);
  print("You provided $nb_name blast name and $nb_file blast files. These number must be equal\n");
  pod2usage(0);
}
#We parse the fasta file and fill the basic information for each feature
my $fasta_handle = Bio::SeqIO->new(
                              -format => 'fasta',
                              -file   => $fasta,
                              );
my %features;

while( my $seq = $fasta_handle->next_seq() )
{
	my $name= $seq->display_name();
	$features{$name} = new Bio::SeqFeature::Generic ( -start => 1, -end => $seq->length(),
                                -strand => 1, -primary => 'mRNA',
                                -source_tag   => 'Abyss_Cap3',
				-seq_id => $prefix."_".$seq->display_name(),
                                -score  => 1,
                                 );
    $features{$name}->add_tag_value("ID",$prefix."_".$name);
    $features{$name}->add_tag_value("Name",$prefix."_".$name);
}
#Parse each XML and start adding Dbxref and annotation fields
for(my $i=0;$i<scalar(@blast_name);$i++)
{
  my $blast_handle = new Bio::SearchIO(-format => 'blastxml', 
                           -file   => $blast_file[$i]);

  while( my $result = $blast_handle->next_result ) 
  {
    my @tag_dbxref;
    my @tag_annotation;
    my $name = $result->query_description;
    while( my $hit = $result->next_hit ) 
    {
      my ($id,$desc) ="unknown ";
      #SP and TrEMBL handling
      if($blast_name[$i] eq "TrEMBL" || $blast_name[$i] eq "SwissProt" )
      {
        my @splitted = split(/\|/,$hit->name);
        $id = $splitted[1];
        $desc = $hit->description;
      }
      #NR and NT handling
      elsif($blast_name[$i] eq "NR" || $blast_name[$i] eq "NT" )
      {
        $id = $hit->accession;
        $desc = $hit->description;
      }
       #MSU and TAIR handling
      elsif($blast_name[$i] eq "MSU" ||$blast_name[$i] eq "TAIR" )
      {
        $id = $hit->name;
        $desc = $hit->description;
      }
      #TAIR handling
      elsif($blast_name[$i] eq "TAIR" )
      {
        my @splitted = split(/\|/,$hit->description);
        $id = $splitted[0];
        $desc = $splitted[2];
      }
     
      $desc=~ tr/,:;/ /;
      
      push(@tag_dbxref, $blast_name[$i].":".$id);
      push(@tag_annotation, $blast_name[$i].":".$desc);
    
    }
    foreach my $dbxref (@tag_dbxref)
    {
      $features{$name}->add_tag_value("Dbxref",$dbxref);
    }
    foreach my $annotations (@tag_annotation)
    {
      $features{$name}->add_tag_value("annotations",$annotations);
    }
  }
                          
}


#Parse Blast2GO and add ontology_term

my $b2g_handle;
if(defined $annot)
{
	if(open($b2g_handle, $annot))
	{
		while(my $line=<$b2g_handle>)
		{
			my @splitted = split(/\t/, $line);
			my $id = $splitted[0];
			my $go = $splitted[1];
			$features{$id}->add_tag_value("Ontology_term",$go);
		}
	}
	else
	{
		print("Cannot open $annot");
		exit(0);
	}
}
#Parse stat and fill cds potential frameshift

my $stat_handle;
if(defined $stat)
{
	if(open($stat_handle, $stat))
	{
		#First parse header
		my %header;
		my $line = <$stat_handle>;
		my @splitted_header = split(/\t/, $line);
		
		for(my $i =0; $i<scalar(@splitted_header); $i++)
		{
			$header{$splitted_header[$i]}=$i;
		}
		
		#Parse other lines to extract information
		while($line = <$stat_handle>)
		{
			my @splitted = split(/\t/, $line);
			my $id = $splitted[$header{"Contig"}];
			my $cds_start = $splitted[$header{"CDS start"}];
			my $cds_stop = $splitted[$header{"CDS end"}];
			my $method = $splitted[$header{"p4e method"}];
			my $shift = "N";
			if((abs($cds_stop - $cds_start)-1)%3 !=0)
			{
				$shift="Y";
			}
			
			$features{$id}->add_tag_value("cds",$cds_start."..".$cds_stop);
			$features{$id}->add_tag_value("p4e method",$method);
			$features{$id}->add_tag_value("frameshift",$shift);
 		}
	}
}

#Write output
my $gff_handle = new Bio::Tools::GFF(-gff_version => 3,
                                     -file => ">".$output);

foreach my $key (keys(%features))
{
  $gff_handle->write_feature($features{$key});
}

$gff_handle->close();

=pod

=head1 DEPENDENCIES

None

=head1 INCOMPATIBILITIES

Fully compatible with any perl version

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=item Gautier SARAH (INRA), gautier.sarah-at-supagro.inra.fr

=back

=head1 VERSION

1.0.0

=head1 DATE

25.06.2014

=head1 LICENSE AND COPYRIGHT

  Copyright 2014 INRA
 
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

