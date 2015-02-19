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

=head1 Shift VCF positions based on an exonerate output file contianing vulgar string


Shift_VCF_positions_based_on_exonerate.pl - Transform VCF transcriptomic coordinate in genomic coordinate based on an exonerate file containing vulgar String

=head1 SYNOPSIS

Shift_VCF_positions_based_on_exonerate.pl -v vcf_file -e exonerate_file -o vcf_file

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
use Data::Dumper;

=pod

=head1 OPTIONS


=head2 Parameters

=over 5

=item B<[-v]> ([a vcf file]): 

The infput VCF file with coordinate to modify 

=item B<[-o]> ([a VCF File]): 

The output VCF file with the coordinate modified

=item B<[-e]> ([a gff file])

A gff file containing the coordinate of the mRNA and its different exons.

=cut

my $version="";

unless (@ARGV)
{
  pod2usage(0);
  exit();
}

#options processing
my ($man, $help, $debug, $input, $output, $exonerate);

# parse options and print usage if there is a syntax error.
#--- see doc at http://perldoc.perl.org/Getopt/Long.html
GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "vcf|v=s" => \$input, 
           "exonerate|e=s"     => \$exonerate,
           "output|o=s"  => \$output,
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

if(!(defined($exonerate) && defined($input) && defined($output)))
{
  pod2usage(0);
  exit();
}

my %transcrit_positions;
my %genome_positions;
my %transcrit_chr;

system("grep vulgar $exonerate >tmp_vulgar");

my $vulgar_handle;

if (open($vulgar_handle, "tmp_vulgar"))
{
  while(<$vulgar_handle>)
  {
    
    my @trans_pos;
    my @genome_pos;
    chomp();
    my @splitted = split(/\s+/);
    my $contig = $splitted[1];
    if(!exists($transcrit_positions{$contig}))
    {
      $transcrit_chr{$contig} = $splitted[5];
      my $genome_strand = $splitted[8];
      my $trans_strand = $splitted[4];
      my $chr_actual = $splitted[6];
      my $trans_actual = $splitted[2];
      
      if($splitted[10] eq 'M')
      {
        push(@trans_pos, $trans_actual);
        push(@genome_pos, $chr_actual);
        
        if($genome_strand eq '+')
        {
          $chr_actual += $splitted[12];
        }
        else
        {
          $chr_actual -= $splitted[12];
        } 
        if($trans_strand eq '+')
        {
          $trans_actual += $splitted[11];
        }
        else
        {
          $trans_actual -= $splitted[11];
        }
        
        
       push(@trans_pos, $trans_actual);
       push(@genome_pos, $chr_actual);
      }
      else
      {
        print "Contig $contig alignment didn't start with a match";
        exit;
      }
      for(my $i=13; $i<scalar(@splitted);$i+=3)
      {
        if($splitted[$i] eq 'M')
        {
          if($genome_strand eq '+')
          {
            $chr_actual += 1;
          }
          else
          {
            $chr_actual -= 1;
          } 
          if($trans_strand eq '+')
          {
            $trans_actual += 1;
          }
          else
          {
            $trans_actual -= 1;
          }
          push(@trans_pos, $trans_actual);
          push(@genome_pos, $chr_actual);
        }
        if($genome_strand eq '+')
        {
          $chr_actual += $splitted[$i+2];
        }
        else
        {
          $chr_actual -= ($splitted[$i+2]);
        }
          
        if($trans_strand eq '+')
        {
          $trans_actual += $splitted[$i+1];
        }
        else
        {
          $trans_actual -=($splitted[$i+1]);
        }
        if($splitted[$i] eq 'M')
        {
          if($trans_strand eq '+')
          {
            push(@trans_pos, --$trans_actual);
          }
          else
          {
            push(@trans_pos, ++$trans_actual);
          }
          if($genome_strand eq '+')
          {
            push(@genome_pos, --$chr_actual);
          }
          else
          {
            push(@genome_pos, ++$chr_actual);
          }
        }
      }
      $transcrit_positions{$contig}=\@trans_pos;
      $genome_positions{$contig} = \@genome_pos;
    }
  }
}
else
{
  print "cannot open tmp_vulgar";
  exit;
}

close($vulgar_handle);

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
        my $contig = $splitted[0];
        my $pos = $splitted[1];
        if(exists($transcrit_positions{$contig}))
        {
          
          my $chr = $transcrit_chr{$contig};
          my @trans_pos = @{$transcrit_positions{$contig}};
          my @chr_pos = @{$genome_positions{$contig}};
         
          $splitted[0] = $chr;
          my $found =0;
          for(my $i = 0;$i<scalar(@trans_pos); $i+=2)
          {
            
            if($pos >= $trans_pos[$i] && $pos <= $trans_pos[$i+1])
            {
              if($chr_pos[$i] < $chr_pos[$i+1])
              {
                $splitted[1] = $chr_pos[$i+1]-($trans_pos[$i+1]-$pos);
                $found++;
              }
              else
              {
                $splitted[1] = $chr_pos[$i+1]+($trans_pos[$i+1]-$pos)+1;
                $splitted[3] =~ tr/acgtACGT/tgcaTGCA/;
                $splitted[4] =~ tr/acgtACGT/tgcaTGCA/;
                $found++;
              }
            }
            elsif($pos <= $trans_pos[$i] && $pos >= $trans_pos[$i+1])
            {
              print $pos;
              print Dumper(\@trans_pos);
              print Dumper(\@chr_pos);
              if($chr_pos[$i] < $chr_pos[$i+1])
              {
                print "non handled case #1 \n $_";
                exit;
              }
              else
              {
                print "non handled case #2 \n $_";
              }
            }
          }
          if($found)
          {
            print $output_handle join("\t",@splitted)."\n";
          }
          else
          {
            print "exists:".exists($transcrit_positions{$contig})."\t".$_."\n";
          }
        }
      }
      else
      {
        print $output_handle $_;
      }
    }
  }
  else
  {
    print "Cannot open $output";
    exit;
  }
}
else
{
  print "Cannot open $input";
  exit;
}
rmtree("tmp_vulgar");




