#!/usr/bin/perl

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

=head1 Shift VCF positions based on gff file


Shift_VCF_positions_based_on_gff.pl - Transform VCF transcriptomic coordinate in genomic coordinate based on a gff file

=head1 SYNOPSIS

Shift_VCF_positions_based_on_gff.pl -v vcf_file -g gff_file -o vcf_file

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
use Bio::Tools::GFF;
use List::Util qw(min max);

=pod

=head1 OPTIONS


=head2 Parameters

=over 5

=item B<[-v]> ([a vcf file]): 

The infput VCF file with coordinate to modify 

=item B<[-o]> ([a VCF File]): 

The output VCF file with the coordinate modified

=item B<[-g]> ([a gff file])

A gff file containing the coordinate of the mRNA and its different exons.

=cut

my $version="";

unless (@ARGV)
{
  pod2usage(0);
  exit();
}

#options processing
my ($man, $help, $debug, $input, $output, $gff);

# parse options and print usage if there is a syntax error.
#--- see doc at http://perldoc.perl.org/Getopt/Long.html
GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "vcf|v=s" => \$input, 
           "gff|g=s"     => \$gff,
           "output|o=s"  => \$output,
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

if(!(defined($gff) && defined($input) && defined($output)))
{
  pod2usage(0);
  exit();
}


my %transcrit_positions;
my %chr_positions;
my %neg_strand_transcrit;
my %chr;

my $gff_in = Bio::Tools::GFF->new(
                                  -file => $gff,
                                  -gff_version => 3
				  );
				  
while(my $feature = $gff_in->next_feature())
{
  if ($feature->primary_tag() eq "CDS" || $feature->primary_tag() =~ /UTR/i || $feature->primary_tag() =~ /exon/i)
  {
   
    my ($transcrit) = $feature->get_tag_values("Parent");
    $transcrit =~ s/\s+$// ;
    $chr{$transcrit}=$feature->seq_id();
    push(@{$chr_positions{$transcrit}}, $feature->start, $feature->end);
    if(!defined($transcrit_positions{$transcrit}))
    {
      push(@{$transcrit_positions{$transcrit}},1,(($feature->end - $feature->start) + 1));
    }
    else
    {
      my $last = $transcrit_positions{$transcrit}[-1];
      push(@{$transcrit_positions{$transcrit}},$last+1,$last+(($feature->end - $feature->start) + 1));
    }
    if($feature->strand == -1)
    {
      $neg_strand_transcrit{$transcrit}++;
    }
  }
  
}

#On retourne les positions du transcrit si sur brin négatif
foreach my $neg_transcript (keys(%neg_strand_transcrit))
{
  my $remaining_length = $transcrit_positions{$neg_transcript}[-1];
  @{$transcrit_positions{$neg_transcript}}=();
  for(my $i=0;$i<@{$chr_positions{$neg_transcript}};$i+=2)
  {
    push(@{$transcrit_positions{$neg_transcript}},$remaining_length,$remaining_length-($chr_positions{$neg_transcript}[$i+1]-$chr_positions{$neg_transcript}[$i]) );
    $remaining_length = $remaining_length-($chr_positions{$neg_transcript}[$i+1]-$chr_positions{$neg_transcript}[$i])-1;
  }
}

#Pour chaque position du VCF, on récupère la position sur le génome et on l'écrit dans la sortie. En inversant les allèles si sur brin négatif
my $input_handle;
my $output_handle;

if (open($input_handle, $input))
{
  if(open($output_handle, ">".$output))
  {
    while(<$input_handle>)
    {
      if($_ !~ /^#/)
      {
        chomp();
        my @splitted = split(/\t/);
        my $trans = $splitted[0];
        my $printed=0;
        if(defined $transcrit_positions{$trans})
        {
          for(my $i=0;$i<@{$transcrit_positions{$trans}};$i+=2)
          {
            if($splitted[1]<=$transcrit_positions{$trans}[$i] && $splitted[1]>=$transcrit_positions{$trans}[$i+1])
            {
              my $chr=$chr{$trans};
              my $chr_pos=$chr_positions{$trans}[$i+1]-($splitted[1]-$transcrit_positions{$trans}[$i+1]);
              $splitted[3] =~ tr/acgtACGT/tgcaTGCA/;
              my $ref = reverse($splitted[3]);
              $splitted[4] =~ tr/acgtACGT/tgcaTGCA/;
              my $var = reverse($splitted[4]);
              $splitted[0]=$chr;
              $splitted[1]=$chr_pos;
              $splitted[3]=$ref;
              $splitted[4]=$var;
              $printed=1;
              print $output_handle join("\t",@splitted)."\n";
            }
            elsif($splitted[1]>=$transcrit_positions{$trans}[$i] && $splitted[1]<=$transcrit_positions{$trans}[$i+1])
            {
              my $chr=$chr{$trans};
              my $chr_pos=$chr_positions{$trans}[$i]+($splitted[1]-$transcrit_positions{$trans}[$i]);
              $splitted[0]=$chr;
              $splitted[1]=$chr_pos;
              $printed=1;
              print $output_handle join("\t",@splitted)."\n";
            }
          }
          if(!$printed)
          {
            print $trans."\t".$_;
            for(my $i=0;$i<@{$transcrit_positions{$trans}};$i++)
            {
              print "\t".$transcrit_positions{$trans}[$i];
            }
            
            print "\n";
          }
        }
      }
    }
  }
}

close $input_handle;
close $output_handle;
