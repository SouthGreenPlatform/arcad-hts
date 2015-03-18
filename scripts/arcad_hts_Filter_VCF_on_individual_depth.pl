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

=head1 Filter VCF


Filter_VCF - Filter a VCF based on individual depth

=head1 SYNOPSIS

Filter_VCF_on_individual_depth -i vcf_file -o vcf_file -ind_prefix individual_prefix -ind_min_depth minimum_depth integer -ind_min_number number -m mask_depth

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Readonly;
use Getopt::Long;
use Pod::Usage;
use Fatal qw/:void open close/;
use File::Path qw (rmtree);


=pod

=head1 OPTIONS


=head2 Parameters

=over 5

=item B<[-i]> ([a vcf file]): 

The input VCF file

=item B<[-o]> ([a VCF File]): 

The output VCF file 

=item B<[-ind_prefix]> ([a string])

A prefix to allowing to groups sample

=item B<[-ind_min_depth]> ([an integer ])

The minimum depth for a sample, to consider it at a position

=item B<[-ind_min_number]> ([an integer])

The minimum number of individuals that are considered at the position to keep it in the output VCF

=item B<[-mask]> ([an integer])

The minimum depth to mask a genotype (cross between a snp position and an individual). All genotypes strictly below this depth will be transformed
by a ./.
Default: 0

=cut

my $version="";

unless (@ARGV)
{
  pod2usage(0);
  exit();
}

#options processing
my ($man, $help, $debug, $input, $output, @prefix, @depth, @number, $mask_depth);

# parse options and print usage if there is a syntax error.
#--- see doc at http://perldoc.perl.org/Getopt/Long.html
GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "ind_prefix|p=s" =>\@prefix,
           "ind_min_depth|d=i" =>\@depth,
           "ind_min_number|n=i" =>\@number,
           "mask|m=i"       =>\$mask_depth,
           "input|conf|i=s" => \$input, 
           "output|o=s"  => \$output,
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

if(@prefix <1)
{
  pod2usage(0);
  exit("You did not precise any sample groups");
}

if(@prefix != @depth)
{
  pod2usage(0);
  exit("You did not precise the same number of sample groups than sample depth");
}

if(@prefix != @number) 
{
  pod2usage(0);
  exit("You did not precise the same number of sample groups than sample number");
}

my $input_handle;
my $output_handle;
my %groups_id;
$|++;
if (open($input_handle, $input))
{
  if(open($output_handle, ">".$output))
  {
    while(<$input_handle>)
    {
      chomp();
      if($_ =~ /^#/)
      {
        if($_ =~ /^#CHROM/)
        {
          my @splitted = split(/\t/);
          foreach my $group (@prefix)
          {
            for(my $i = 9;$i<@splitted;$i++)
            {
              if($splitted[$i] =~ /^$group/)
              {
                push(@{$groups_id{$group}},$i);
              }
              
            }
          }
        }
        print $output_handle $_."\n";
      }
      else
      {
        my $i=0;
        my $ok =1;
        my @splitted = split(/\t/);
        #if($splitted[6] eq "PASS")
        #{
          foreach my $group (@prefix)
          {
            my $num =0;
            
            foreach my $id (@{$groups_id{$group}})
            {
              if($splitted[$id] ne "./.")
              {
                my @split_id = split(':', $splitted[$id]);
                if($split_id[2] >= $depth[$i])
                {
                  $num++;
                }
                if($split_id[2]<$mask_depth)
                {
                  $splitted[$id] = "./.";
                }
              }
            }
            if($num<$number[$i])
            {
              $ok=0;
            }
            $i++;
          }
          if($ok)
          {
            print $output_handle join("\t", @splitted)."\n";
          }
        #}
        
      }
    }
  }
}
else
{
	print "Cannot open $input";
	exit();
}
