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

=head1 VCF to Fasta

Reconstruct a fasta from a VCF file. Reconstruct a fasta file for each sample specified.

=head1 SYNOPSIS

VCF_to_Fasta.pl -v vcf_file -o fasta_file -r fasta_reference -c cov_threshold -b bam_file -ind_prefix individual_prefix -ind_min_depth integer -L intervals_file 

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

use lib '/NAS/arcad_data/Softs/tags';
use Modules::Config::Softwares;
use Modules::Files::Bam;

=pod

=head1 OPTIONS


=head2 Parameters

=over 5

=item B<[-v]> ([a vcf file]): 

The input VCF file containing SNPs. Only SNP with flag passed will be used

=item B<[-b]> ([a bam file]): 

The input BAM file used for the SNP Calling. Only SNP with flag PASS will be used

=item B<[-c]> ([an integer): 

Coverage trheshold. Positions having coverage below this threshold in a sample will be replace by a N in this sample's fasta

=item B<[-o]> ([a directory]): 

The output directory. This directory will contain the output fasta file. Fasta file will be named as follow <sample>.fasta

=item B<[-r]> ([a fasta File]): 

The reference fasta file. The reference file used for the mapping and the SNP_Calling

=item B<[-ind_name]> ([a string])

Sample name for which a fasta file will be created. Can be specified more than once.
If this is not set, then all samples will be extracted.

=item B<[-L]> ([a interval file]): 

An interval file in order to restrain the output to some positions. Each line must be tabulated like this: chr<tab>start<tab>end<tab>name (the name is optional, if not provided, name will be number)

=cut

my $version="";

unless (@ARGV)
{
  pod2usage(0);
}

    my %IUPAC = (
        A    => 'A',
        T    => 'T',
        U    => 'U',
        C    => 'C',
        G    => 'G',
        AC   => 'M',
        CA   => 'M',
        AG   => 'R',
        GA   => 'R',
        AT   => 'W',
        TA   => 'W',
        CG   => 'S',
        GC   => 'S',
        CT   => 'Y',
        TC   => 'Y',
        GT   => 'K',
        TG   => 'K'
    );

my $SAMTOOLS_PATH=&$Softwares::SAMTOOLS_PATH or confess("$!");
my $INTERSECTBED_PATH=&$Softwares::INTERSECTBED_PATH or confess("$!");

#options processing
my ($man, $help, $debug, $vcf, $bam, $output, $ref, @prefix, $depth, $intervals);
$depth = 10;
$debug = 0;
# parse options and print usage if there is a syntax error.
#--- see doc at http://perldoc.perl.org/Getopt/Long.html
GetOptions("help|?"     => \$help,
           "man"        => \$man,
           "debug:i"    => \$debug,
           "ind_prefix|p=s" =>\@prefix,
           "coverage|c=i" =>\$depth,
           "vcf|v=s" => \$vcf, 
           "reference|r=s" => \$ref,
           "intervals|L=s" => \$intervals, 
           "bam|b=s"  => \$bam,
           "output|o=s"  => \$output,
) or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}

if(!defined($ref))
{
  print("You did not specify a reference\n");
  pod2usage(0);
}
if(! -e $ref) 
{
  print("Cannot find $ref");
  pod2usage(0);
}

if(! -e $output)
{
  mkdir $output;
}

# INDEXING input bam file
my @t = ($bam);
my $o_bam = Bam->new( \@t );
$o_bam->indexBam;

my $vcf_handle;
my $output_handle;
my %groups_id;
$|++;
my $tmp = "tmp_".$$;
mkdir $tmp;

#If no prefix were indicated, all samples from the VCF will be extracted
if(@prefix <1)
{
  if (open($vcf_handle, $vcf))
  {
    while(<$vcf_handle>)
    {
      if($_ =~ /^#/)
	  {
        if($_ =~ /^#CHROM/)
	    {
		  chomp();
	      my @splitted = split(/\t/);  
		  for(my $i = 9;$i<@splitted;$i++)
		  {
		    my $group = $splitted[$i];
		    push(@prefix, $group);
		    push(@{$groups_id{$group}},$i);  
		  }
        }
	  }
    }
  }
  else
  {
    print "Cannot open $vcf";
    exit();
  }
}

#Get sample positions in the VCF file
elsif (open($vcf_handle, $vcf))
{
  while(<$vcf_handle>)
  {
    if($_ =~ /^#/)
      {
        if($_ =~ /^#CHROM/)
        {
          my @splitted = split(/\t/);
          foreach my $group (@prefix)
          {
            for(my $i = 9;$i<@splitted;$i++)
            {
              if($splitted[$i] =~ /^$group$/)
              {
                push(@{$groups_id{$group}},$i);
              }
            }
          }
        }
      }
  }
}
else
{
  print "Cannot open $vcf";
  exit();
}
close($vcf_handle);


#Create a vcf file focused on the given intervals

if(defined($intervals))
{
	system("$INTERSECTBED_PATH -a $vcf -b $intervals >$tmp/intervals.vcf");
}


#Create the fasta for each given intervals

if(! -e "$ref.fai")
{
  system("$SAMTOOLS_PATH faidx $ref");
}

system("sed 's/\t/:/' $intervals | sed 's/\t/-/' |sed 's/\t.*//' >$tmp/intervals");

mkdir "$tmp/fasta";

my $intervals_handle;
my $numerical_id=0;
if (open($intervals_handle, "$intervals"))
{
  while(<$intervals_handle>)
  {
     chomp();
     my @splitted = split('\t');
     my $interval = $splitted[0].":".$splitted[1]."-".$splitted[2];

     my $fasta_file = @splitted == 3 ? $numerical_id++ : $splitted[3];

	system("$SAMTOOLS_PATH faidx $ref '$interval' >$tmp/fasta/$fasta_file.fasta && sed -i 's/>.*/>'$fasta_file'/' $tmp/fasta/$fasta_file.fasta");
  }
}
else
{
  print "Cannot open $intervals";
  exit();
}
close($intervals_handle);

$numerical_id=0;
#For each sample

foreach my $group (@prefix)
{
  my $fasta_handle;
  open($fasta_handle, ">$output/$group.fasta");
  
#  For each interval
  if (open($intervals_handle, "$intervals"))
  {
    while(<$intervals_handle>)
    {
      chomp();
      my @splitted = split('\t');
      my $interval = $splitted[0].":".$splitted[1]."-".$splitted[2];

      my $name = @splitted == 3 ? $numerical_id++ : $splitted[3];
      #Recover original fasta
      my $seqio_object = Bio::SeqIO->new(-file => "$tmp/fasta/$name.fasta");
      my $seq = $seqio_object->next_seq;
      $seq->display_id($name);
      my $nucleotide = $seq->seq;
      #getSNPs and modify the sequence
      my $vcf_handle;
      if(open($vcf_handle, "$tmp/intervals.vcf"))
      {
        my $id = @{$groups_id{$group}}[0];
        while(<$vcf_handle>)
        {
        chomp();
        my @vcf_col = split("\t");

        if($vcf_col[1] <= $splitted[2] && $vcf_col[1] >= $splitted[1] && ($vcf_col[6] eq "PASS" || $vcf_col[6] eq "DP_FILTER"))
        {
          if(length($vcf_col[3]) == 1 && length($vcf_col[4]) == 1)
          {
            if($vcf_col[$id] ne "./.")
            {
              my @field = split(":", $vcf_col[$id]);
              my $pos = $vcf_col[1] - $splitted[1] ;
              if($field[0] eq "0/1" || $field[0] eq "0|1")
              {
                substr($nucleotide, $pos, 1) = $IUPAC{$vcf_col[3].$vcf_col[4]};
              }
              elsif($field[0] eq "1/1" || $field[0] eq "1|1")
              {
                substr($nucleotide, $pos, 1) = $vcf_col[4];
              }
            }
          }
        }
      }
     }
     else
     {
       print "Cannot open $tmp/intervals.vcf";
       exit();
     }
     close($vcf_handle);
   
   #Create BAM
	 if($debug >0)
	 {
		print "Launched:\n"."$SAMTOOLS_PATH view -r $group -h -b  $bam '$interval'>$tmp/tmp.bam && $SAMTOOLS_PATH index $tmp/tmp.bam\n";
	 }
     system("$SAMTOOLS_PATH view -r $group -h -b  $bam '$interval'>$tmp/tmp.bam && $SAMTOOLS_PATH index $tmp/tmp.bam");
   #Create depth
     system("$SAMTOOLS_PATH depth -r '$interval' $tmp/tmp.bam >$tmp/tmp.dp");
   #If depth is empty, we look for the ID in the SM fields of RG
     if(-z "$tmp/tmp.dp")
     {
       system("$SAMTOOLS_PATH view -H  $bam | sed -n '/^".'@RG'."/ p' >$tmp/RG.tmp");
       my $rg_handle;
       my $rg_out;
       my @rg_to_keep;
       if(open($rg_handle, "$tmp/RG.tmp") && open($rg_out, ">$tmp/rg_list"))
	   {
		 while(<$rg_handle>)
		 {
	       chomp();
	       if($_ =~ /ID:([^\s]+).*SM:$group/)
	       {
		       print $rg_out $1."\n";
	       }
	     }
	   }
     }
	 system("$SAMTOOLS_PATH view -R $tmp/rg_list -h -b  $bam '$interval'>$tmp/tmp.bam && $SAMTOOLS_PATH index $tmp/tmp.bam");
	 system("$SAMTOOLS_PATH depth -r '$interval' $tmp/tmp.bam >$tmp/tmp.dp");
   #Replace by N where depth is below threshold
     my $depth_handle;
   #If the depth file is empty, then we replace all the sequence by N's
     if(-e "$tmp/tmp.dp" && -z "$tmp/tmp.dp")
     {
		 $nucleotide =~ tr/ACGTMRWSYKacgtmrwsyk/N/ ;
	 }
     elsif(open($depth_handle,"$tmp/tmp.dp" ))
     {
        my $line = <$depth_handle>;
        my @dp_field = split("\t", $line);
        my $pos = $dp_field[1] - $splitted[1] ;
       #if the first line of the depth file does not correpond to the
       #beginning of the sequence, then replace all the beginning by N's
       if( $pos != 0)
        {
          for(my $i=0; $i<$pos;$i++)
          {
            substr($nucleotide,$i,1) = "N";
          }
          if($dp_field[2] < $depth)
          {
            substr($nucleotide,$pos,1) = "N";
          }
        }
        my $old = $pos;
        while($line = <$depth_handle>)
        {
          @dp_field = split("\t", $line);
          $pos = $dp_field[1] - $splitted[1] ;
          if(($pos - $old) > 1)
          {
            for(my $i=$old+1; $i<$pos;$i++)
            {
              substr($nucleotide,$i,1) = "N";
            }
            if($dp_field[2] < $depth)
            {
              substr($nucleotide,$pos,1) = "N";
            }
          }
          elsif(($pos - $old) == 1)
          {
            if($dp_field[2] < $depth)
            {
              substr($nucleotide,$pos,1) = "N";
            }
          }
          else
          {
            print "Problem in parsing depth file $pos $old ".$dp_field[1];
            exit();
          }
          $old = $pos;
        }
        close($depth_handle);
        for(my $i=$pos+1; $i<(length($nucleotide));$i++)
          {
            substr($nucleotide,$i,1) = "N";
          }
      }
      else
      {
        print "Cannot open $tmp/tmp.dp";
        exit();
      }
      
      print $fasta_handle ">".$name."\n".$nucleotide."\n";
    }
  }
  close($intervals_handle);
  close($fasta_handle);
}


=head1 DEPENDENCIES

Samtools, BedTools (intersectBed)

=head1 INCOMPATIBILITIES

Fully compatible with any perl version

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 2

=item Gautier SARAH (INRA), gautier.sarah-at-supagro.inra.fr

=back

=head1 VERSION

1.0.0

=head1 DATE

23/01/2014

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
