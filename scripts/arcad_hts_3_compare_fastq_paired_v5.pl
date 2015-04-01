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

#$ -q arcad.q
#$ -cwd
#$ -V
#$ -b yes
#$ -N Compare_fastq_paired


#This script will compare two FASTQ files and check their homogeneity in term of forward-reverse.
#It is based on their name
#The non-paired sequences will be printed in a third file, and can be treated latter as single sequences.
#The files contain pairs but are not truely paired, ie the first sequence of the 1st file did not correspond to the first on the 2nd file...
#It will work for Fastq with one line for sequence, one line for qual ONLY

=pod

=head1 NAME

compare_fastq_paired_v5 - Resynchronized Fastq paired files desynchonized by cleaning

=head1 DESCRIPTION

Resynchronized Fastq paired files desynchonized by cleaning

=head1 SYNOPSIS / USAGE

compare_fastq_paired_v5 -f R1_file -r R2_file [-of R1_output] [-or R2_output] [-os unpaired_output]  

=cut

use strict;
use warnings;
use Getopt::Long;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Time::HiRes qw( time );
use Pod::Usage;

#print the line when they are in buffer
$|=1;

my $start_time =time();

pod2usage(0) unless (@ARGV);

=head1 OPTIONS

=head2 required parameters

=over 1

=item B<-f > (R1 Fastq file)

Fastq file to proccess (can be fastq.gz)

=item B<-r > (R2 Fastq file)

Fastq file to proccess (can be fastq.gz)

=back

=head2 optional arguments

=over 1

=item B<-of > (R1 output Fastq file)

Fastq file to proccess (can be fastq.gz)

=item B<-or > (R2 output Fastq file)

Fastq file to proccess (can be fastq.gz)

=item B<-os > (unpaired Fastq file)

Fastq file to proccess (can be fastq.gz)

=item B<-? | --man | --help>

Print this help

=back

=cut

my $version="";	


my ($forward,$reverse,$output_forward,$output_reverse,$output_single,$help, $man);

GetOptions("help|?|h" => \$help,
		"man"        => \$man,
		"f|forward=s"=>\$forward,
		"r|reverse=s"=>\$reverse,
		"of|output_forward=s"=>\$output_forward,
		"or|output_reverse=s"=>\$output_reverse,
		"os|output_single=s"=>\$output_single
		)or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}


my ($forwardname,$reversename, $singlename);
#Creating output names:
if(!defined($output_forward))
{
  my @forwardlist=split/\//,$forward;
  $forwardname=$forwardlist[-1];
  $forwardname=~s/\.\w{1,}(\.gz)?$//;
  $forwardname.="_paired.fastq";
  if(defined($1))
  {
	  $forwardname.=$1;
  }
  
}
else
{
  $forwardname=$output_forward;
}

if(!defined($output_reverse))
{
  my @reverselist=split/\//,$reverse;
  $reversename=$reverselist[-1];
  $reversename=~s/\.\w{1,}(\.gz)?$//;
  $reversename.="_paired.fastq";
  if(defined($1))
  {
	  $reversename.=$1;
  }
}
else
{
  $reversename=$output_reverse;
}

if(!defined($output_single))
{
  my @forwardlist=split/\//,$forwardname;
  $singlename=$forwardlist[-1];
  $singlename=~s/\.\w{1,}(\.gz)?$//;
  $singlename.="_single.fastq";
  if(defined($1))
  {
	  $singlename.=$1;
  }
}
else
{
  $singlename=$output_single;
}

print("Creating hash \n");
#Gathering name and modifying them to be able to compare.
my $comheadline;
if($forward =~ m/\.gz$/)
{
	$comheadline="zcat $forward |head -n 1 ";
}
else
{
	$comheadline="head -n 1 $forward";
}
my $headline=`$comheadline`;
chomp $headline;
$headline=substr($headline,0,4);

my $comid;
if($forward =~ m/\.gz$/)
{
	$comid="zcat $forward | sed -n '/$headline/ p' ";

}
else
{
	$comid="sed -n '/$headline/ p' $forward";

}
my $list_id_for= `$comid`;
chomp $list_id_for;
my @id_for=split/\n/,$list_id_for;
#print $list_id_for ;
#free some perl memory
undef($list_id_for);

my %idhash;
foreach my $name (@id_for)
{
  $name=~s/\/\d$//;
  $name=~s/\s\d:.+$//;
  $idhash{$name}++;
  
}

my $duration=time() - $start_time;
print("Hash created. Duration:". getDuration($duration)." \n");

#free some perl memory
undef(@id_for);
print(system("ps aux | grep $$"));

#Reading and outputing
open(my $single_out_handle,">$singlename") or die("\nCannot create $singlename file: $!\n");
if($singlename =~ m/\.gz$/)
{
	$single_out_handle = new IO::Compress::Gzip $single_out_handle or die "IO::Compress::Gzip failed: $GzipError\n";
	$single_out_handle->autoflush(1);
}
my $reverse_handle;
if($reverse =~ m/\.gz$/)
{
	$reverse_handle = new IO::Uncompress::Gunzip $reverse or die "gunzip failed: $GunzipError\n";
}
else
{
	open($reverse_handle,"$reverse") or die("\nCannot read $reverse file: $!\n");

}
my $reverse_out_handle;
if($reversename =~ m/\.gz$/)
{
	$reverse_out_handle = new IO::Compress::Gzip $reversename or die "IO::Compress::Gzip failed: $GzipError\n";
	$reverse_out_handle->autoflush(1);
}
else
{
	open($reverse_out_handle,">$reversename") or die("\nCannot create $reversename file: $!\n");

}
print "\nPrinting rev...\n";
my %id_both;
while (<$reverse_handle>) #Writing forward correct
    {
    my $in=$_;
    chomp $in;
    my $idhere=$in;
    $idhere=~s/\/\d$//;
    $idhere=~s/\s\d:.+$//;
    $in.="\n".<$reverse_handle>.<$reverse_handle>.<$reverse_handle>;#Pick up the 4 lines
    if (exists($idhash{$idhere})) #They are in the common list
        {
        print $reverse_out_handle $in;
        delete($idhash{$idhere});
        $id_both{$idhere}++;        
        }
    else #They are not in the common list
        {
        print $single_out_handle $in;
        }
    }
close $reverse_handle;
close $reverse_out_handle;
$duration=time() - $start_time;
print("Reverse written. Duration: ". getDuration($duration)." \n");
print(system("ps aux | grep $$"));

my $forward_handle;
if($forward =~ m/\.gz$/)
{
	$forward_handle = new IO::Uncompress::Gunzip $forward or die "gunzip failed: $GunzipError\n";
}
else
{
	open($forward_handle,"$forward") or die("\nCannot read $forward file: $!\n");

}
my $forward_out_handle;
if($forwardname =~ m/\.gz$/)
{
	$forward_out_handle = new IO::Compress::Gzip $forwardname or die "IO::Compress::Gzip failed: $GzipError\n";
	$forward_out_handle->autoflush(1);
}
else
{
	open($forward_out_handle,">$forwardname") or die("\nCannot create $forwardname file: $!\n");

}

print "\nPrinting for...\n";
while (<$forward_handle>) #Writing reverse correct
    {
    my $in=$_;
    chomp $in;
    my $idhere=$in;
    $idhere=~s/\/\d$//;
    $idhere=~s/\s\d:.+$//;
    $in.="\n".<$forward_handle>.<$forward_handle>.<$forward_handle>; #Pick up the 4 lines
    if ($id_both{$idhere}) #They are in the common list
        {
        print $forward_out_handle $in;
        }
    else #They are not in the common list
        {
        if(!exists($idhash{$idhere}))
          {
           print "Check $idhere\n";
          }
        print $single_out_handle $in;
        }
    }
close $forward_handle;
close $forward_out_handle;
close $single_out_handle;

$duration=time() - $start_time;
print("Forward written. Duration: ". getDuration($duration)." \n");
print(system("ps aux | grep $$"));


############################
#
#  SUB
#
############################


sub getDuration
{
    my $duration = shift;
    if (60 > $duration)
    {
        # less than a minute
        return sprintf('%4.2f second(s)', $duration);
    }
    elsif (3600 > $duration)
    {
        # less than an hour
        return sprintf('%d:%04.2f', int($duration/60), ($duration%60));
    }
    elsif (86400 > $duration)
    {
        # less than a day
        return sprintf('%d:%02d:%04.2f', int($duration/3600), int($duration/60)%60, ($duration%60));
    }
    else
    {
        # more than a day
        return sprintf('%d days and %d:%02d:%04.2f', int($duration/86400), int($duration/3600)%24, int($duration/60)%60, ($duration%60));
    }
}


sub modifname { # Remove the 1 or 2 (or any text separated by a space from the name) at the end of the name
    my @list=@_;
    my @out;
    foreach my $name (@list)
        {
        $name=~ s/\/\d$//; # Remove /1 or /2 for Illumina1.3 or 1.5
        $name =~ s/\s\d:.+$//;# Remove 1:..:...: or anything after a space, for Illumina1.8+
        push @out, $name;
        }
    return @out;
    }
    
exit(0);
