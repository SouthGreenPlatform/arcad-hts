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

=head1 Classify blast output

Classify_Blast-  Classify Blast Queries in differents categories

=head1 SYNOPSIS

Classify_Blast.pl -i input -o directory -identity float -merge_hsp -blastx

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Getopt::Long;
use Pod::Usage;
use Fatal qw/:void open close/;
use Graph;
use List::Util qw[min max];

=pod

=head1 OPTIONS


=head2 Parameters

=over 3

=item B<[-i]> ([a blast output]): 

This file must be tabulated. 
First column must contain the query name
Second column must contain the query length
Third column must contain the target name
Fourth column must contain target length
Fifth column must contain the alignment length
Sixth column must contain the identity ratio of the alignment
Seventh and eighth must respectively contain the position of the start and the end of the alignment on the query
Nineth and tenth must respectively contain the position of the start and the end of the alignment on the target
Eleventh must be the evalue

=item B<[-d]> ([a directory]): 

The output directory, where to write the output files.
Seven output files will be produced.
A summary of the repartition of the queries in the different categories.
For each category, a file containing the list of queries in it.
Categories are, Full, Partial, Fragment, Alleles, Chimere, Multicopy

=item B<[-m/nom]> ([a flag])

Determine if HSPs should be merged. If yes, HSP will be sorted on their alignment length.
All HSP on the opposite strain of the firs one will be ignored.
Default: HSPs are merged

=item B<[-blastx]> (a flag)

Determine if the blast was a blastX in order to correct the alignment length, and the coordinate.
Default: blastx is off

=item B<[-id]> ([a float])

Identity required (in pourcentage) to consider a Hit.
Default:80

=cut

my $version="";

unless (@ARGV)
{
  pod2usage(0);
  exit();
}

#options processing
my ($man, $help, $debug, $input, $directory ,$merge_hsp, $identity , $blastx);
$merge_hsp=1;
$identity=80;
# parse options and print usage if there is a syntax error.
#--- see doc at http://perldoc.perl.org/Getopt/Long.html
GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "input|i:s"  => \$input,
           "id:f"       => \$identity,
           "blastx|x!"  => \$blastx,
           "directory|d:s" => \$directory,
           "merge_hsp|m!"  => \$merge_hsp
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

#Creation of output directory if necessary
if(-e $directory)
{
  if( !(-d $directory))
  {
    warn "$directory already exists but it's not a directory";
    exit();
  }
}
else
{
  mkdir $directory;
}

#Creation of a hash containing the characteristic of each Hit (Merging HSP if needed)

my %hits;

#Sort the blast output by descending alignment length
if($merge_hsp)
{
  system("sort -k5,5nr $input >${input}_sorted");
}


#Parsing of the blast output

my $blast_handle;

if (open($blast_handle, "${input}_sorted"))
{
	while(<$blast_handle>)
	{
	  chomp();
	  if($_!~ /^\s*$/)
	  {
	    my @splitted = split(/\t/);
	    if($splitted[5]>=$identity)
	    {
	      if($blastx)
	      {
	        $splitted[3]*=3;
	        $splitted[4]*=3;
	        $splitted[8]*=3;
	        $splitted[9]*=3;
	      }
	      my $id = $splitted[0]."->".$splitted[2];
	      if(exists($hits{$id}) && $merge_hsp)
	      {
	        my @previous = split(/\t/,$hits{$id});
	        
          $previous[4]+= $splitted[4];
          if($splitted[8]>$splitted[9] && $previous[8]>$previous[9])
          {
            $previous[8] = max($previous[8],$splitted[8]);
            $previous[9] = min($previous[9],$splitted[9]);
            
            $previous[6] = min($previous[6],$splitted[6]);
            $previous[7] = max($previous[7],$splitted[7]);
          }
          elsif($splitted[8]<$splitted[9] && $previous[8]<$previous[9])
          {
            $previous[8] = min($previous[8],$splitted[8]);
            $previous[9] = max($previous[9],$splitted[9]);
            
            $previous[6] = min($previous[6],$splitted[6]);
            $previous[7] = max($previous[7],$splitted[7]);
          }
          else
          {
            print "$id seems to have HSP on both strand: Reversed HSP not taken into account\n";
          }
          $hits{$id}=join("\t",@previous);
	        
	        
	      }
	      elsif(!exists($hits{$id}))
	      {
	        $hits{$id}=$_;
	      }
	    }
	  }
	  
	}
	
}
close($blast_handle);
#Creation of the grap containing the connections between query and target
my $connections = Graph->new();
foreach my $hit (keys(%hits))
{
  my @ids = split(/->/,$hit);
  my @splitted = split(/\t/, $hits{$hit});
  if($splitted[4]>= min($splitted[1],$splitted[3])*0.5)
  {
    $connections->add_edge($ids[0],$ids[1]);
  }
  
}

#Parsing of the graph

my @subgraphs = $connections->weakly_connected_components();
my (@full, @partial, @fragment, @allele, @chimere, @multicopy);

foreach my $connected_component (@subgraphs)
{
  my @queries;
  my @targets;
  foreach my $vertex (@{$connected_component})
  {
    if($connections->is_source_vertex($vertex))
    {
      push(@queries,$vertex);
    }
    elsif($connections->is_sink_vertex($vertex))
    {
       push(@targets,$vertex);
    }
    else
    {
      print "Check $vertex, as it seems to be a query AND a target\n";
    }
  }
  
  #Categorizing each subgraph
  if(scalar(@targets) == 1 && scalar(@queries) == 1)
  {
    my @hit = split(/\t/,$hits{$queries[0]."->".$targets[0]});
    
    if($hit[4]>=$hit[3]*0.9)
    {
      push(@full,@queries);
    }
    else
    {
      push(@partial,@queries);
    }
  }
  elsif(scalar(@targets) ==1)
  {
    #We compare blast position of the queries on the targets to check if they overlap above 50% of the smallest query length
    my %alleles;
    for(my $i=0;$i<scalar(@queries)-1;$i++)
    {
      for(my $j=$i+1;$j<scalar(@queries);$j++)
      {
        my $id_1 = $queries[$i]."->".$targets[0];
        my @splitted_1 = split(/\t/, $hits{$id_1});
        my $id_2 = $queries[$j]."->".$targets[0];
        my @splitted_2 = split(/\t/, $hits{$id_2});
        
        my $start_1 = min($splitted_1[8],$splitted_1[9]);
        my $end_1   = max($splitted_1[8],$splitted_1[9]);
        my $start_2 = min($splitted_2[8],$splitted_2[9]);
        my $end_2   = max($splitted_2[8],$splitted_2[9]);
        
        
        if(min($splitted_1[1], $splitted_2[1])*0.5 <= min($end_1,$end_2)-max($start_1,$start_2))
        {
          $alleles{$queries[$i]}++;
          $alleles{$queries[$j]}++;
        }
      }
    }
    
    foreach my $query (@queries)
    {
      if(exists($alleles{$query}))
      {
        push(@allele,$query);
      }
      else
      {
        push(@fragment,$query);
      }
      
    }
  }
  elsif(scalar(@queries) ==1)
  {
     #We compare blast position of the targets on the queries to check if they overlap above 50% of the smallest query length
    my %multicopy;
    for(my $i=0;$i<scalar(@targets)-1;$i++)
    {
      for(my $j=$i+1;$j<scalar(@targets);$j++)
      {
        my $id_1 = $queries[0]."->".$targets[$i];
        my @splitted_1 = split(/\t/, $hits{$id_1});
        my $id_2 = $queries[0]."->".$targets[$j];
        my @splitted_2 = split(/\t/, $hits{$id_2});
        
        my $start_1 = min($splitted_1[8],$splitted_1[9]);
        my $end_1   = max($splitted_1[8],$splitted_1[9]);
        my $start_2 = min($splitted_2[8],$splitted_2[9]);
        my $end_2   = max($splitted_2[8],$splitted_2[9]);
        
        
        if(min($splitted_1[1], $splitted_2[1])*0.5 <= min($end_1,$end_2)-max($start_1,$start_2))
        {
          $multicopy{$targets[$i]}++;
          $multicopy{$targets[$j]}++;
        }
      }
    }
    
    if(scalar(keys(%multicopy)>0))
    {
      push(@multicopy,$queries[0]);
    }
    else
    {
      push(@chimere, $queries[0]);
    }
  }
  else
  {
    if(scalar(@targets) == scalar(@queries))
    {
      my %already_seen;
      my @potential_full;
      my @potential_partial;
      foreach my $query (@queries)
      {
        my @hit = $connections->successors($query);
        my @best_hit;
        foreach my $result (@hit)
        {
          my $id = $query."->".$result;
          my @splitted = split(/\t/, $hits{$id});
          if(scalar(@best_hit)==0)
          {
            @best_hit=@splitted;
          }
          else
          {
            if($splitted[10]>$best_hit[10])
            {
              @best_hit=@splitted;
            }
            elsif($splitted[10] == $best_hit[10] && $splitted[4]>$best_hit[4])
            {
              @best_hit=@splitted;
            }
          }
          
        }
        if($best_hit[4]>=$best_hit[3]*0.9)
        {
          push(@potential_full,$query);
        }
        else
        {
          push(@potential_partial,$query);
        }
        $already_seen{$best_hit[2]}++;
      }
      if(scalar(keys(%already_seen)) == scalar(@targets))
      {
        push(@full, @potential_full);
        push(@partial, @potential_partial);
      }
      else
      {
        push(@multicopy,@queries);
      }
    }
    else
    {
      push(@multicopy,@queries);
    }
    
  }

}
  
print "Full ".scalar(@full)."\n";
print "Partial ".scalar(@partial)."\n";
print "Fragment ".scalar(@fragment)."\n";
print "Allele ".scalar(@allele)."\n";
print "Chimere ".scalar(@chimere)."\n";
print "Multicopy ".scalar(@multicopy)."\n";

my $handle;
open($handle,">$directory/full.lst");
print $handle join("\n", @full);
close($handle);

open($handle,">$directory/partial.lst");
print $handle join("\n", @partial);
close($handle);

open($handle,">$directory/fragment.lst");
print $handle join("\n", @fragment);
close($handle);

open($handle,">$directory/allele.lst");
print $handle join("\n", @allele);
close($handle);

open($handle,">$directory/chimere.lst");
print $handle join("\n", @chimere);
close($handle);

open($handle,">$directory/multicopy.lst");
print $handle join("\n", @multicopy);
close($handle);


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








