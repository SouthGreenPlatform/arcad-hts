#! /bin/env perl

#  
#  Copyright 2015 INRA-CIRAD-IRD
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
#$ -N synchronized_paired_fastq

=pod

=head1 NAME

compare_fastq_paired_v5 - Resynchronized Fastq paired files desynchonized by cleaning

=head1 DESCRIPTION

Resynchronized Fastq paired files desynchonized by cleaning. 
The files need to have the reads on four lines (name,nucleotidic sequence,+, quality sequence)
The files have to be in the same order as their original file, with just missing or trimmed reads


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

#Opening input files
my $reverse_handle;
if($reverse =~ m/\.gz$/)
{
	$reverse_handle = new IO::Uncompress::Gunzip $reverse or die "gunzip failed: $GunzipError\n";
}
else
{
	open($reverse_handle,"$reverse") or die("\nCannot read $reverse file: $!\n");

}

my $forward_handle;
if($forward =~ m/\.gz$/)
{
	$forward_handle = new IO::Uncompress::Gunzip $forward or die "gunzip failed: $GunzipError\n";
}
else
{
	open($forward_handle,"$forward") or die("\nCannot read $forward file: $!\n");

}

#Opening output files
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

open(my $single_out_handle,">$singlename") or die("\nCannot create $singlename file: $!\n");
if($singlename =~ m/\.gz$/)
{
	$single_out_handle = new IO::Compress::Gzip $single_out_handle or die "IO::Compress::Gzip failed: $GzipError\n";
	$single_out_handle->autoflush(1);
}

my $end_of_file=0; #This variable will be set to 1 if we reach the end of the forward input and to 2 if we reach the end of the reverse input
my %forward_name_list;
my %reverse_name_list;

my @forward_read;
my @reverse_read;

my $deleted_forw = 0;
my $deleted_rev = 0;

my $inserted_forw = 0;
my $inserted_rev = 0;


while(!$end_of_file)
{
	if(my $id_forw = <$forward_handle>)
	{ 
		my $read_forw .= $id_forw;
		chomp($id_forw);
		$id_forw=~s/\/\d$//;
		$id_forw=~s/\s\d.+$//;
		$read_forw .= <$forward_handle>.<$forward_handle>.<$forward_handle>;
		
		if(my $id_rev = <$reverse_handle>)
		{
			my $read_rev .= $id_rev;
			chomp($id_rev);
			$id_rev=~s/\/\d$//;
			$id_rev=~s/\s\d.+$//;
			$read_rev .= <$reverse_handle>.<$reverse_handle>.<$reverse_handle>;
			
			if($id_rev eq $id_forw)
			{
				print $forward_out_handle $read_forw;
				print $reverse_out_handle $read_rev;
			}
			else
			{
				#First we check the existence of the reads ids in the opposite name_list array
				my $pos;
				my $found = 0;
				
				if(exists $forward_name_list{$id_rev})
				{
					$pos = $forward_name_list{$id_rev};
					#we flush all the previous reads in the single
					for(my $i=0; $i<$pos-$deleted_forw;$deleted_forw++)
					{
						my $to_print = shift(@forward_read);
						print $single_out_handle $to_print;
						my @splitted = split("\n",$to_print);
						$to_print = $splitted[0];
						$to_print =~s/\/\d.*$//;
						$to_print =~s/\s\d.+$//;
						delete($forward_name_list{$to_print});
					}
					
					for(my $i=scalar(@reverse_read); $i>0;$i--)
					{
						$deleted_rev++;
						my $to_print = shift(@reverse_read);
						print $single_out_handle $to_print;
						my @splitted = split("\n",$to_print);
						my $id = $splitted[0];
						$id =~s/\/\d.*$//;
						$id =~s/\s\d.+$//;
						delete($reverse_name_list{$id});
					}
					
					
					#We flush the matching reads in their respective file
					print $reverse_out_handle $read_rev;
					my $to_print = shift(@forward_read);
					$deleted_forw++;
					print $forward_out_handle $to_print;
					my @splitted = split("\n",$to_print);
					$to_print = $splitted[0];
					$to_print =~s/\/\d.*$//;
					$to_print =~s/\s\d.+$//;
					delete($forward_name_list{$to_print});
					
					
					#We add the forward read to the arrays
					$forward_name_list{$id_forw} = $inserted_forw;
					$inserted_forw++;
					push(@forward_read, $read_forw);
					
					$found = 1;
				}
				
				if(exists $reverse_name_list{$id_forw})
				{
					$pos = $reverse_name_list{$id_forw};
					#we flush all the previous reads in the single
					for(my $i=0; $i<$pos-$deleted_rev;$deleted_rev++)
					{
						my $to_print = shift(@reverse_read);
						print $single_out_handle $to_print;
						my @splitted = split("\n",$to_print);
						$to_print = $splitted[0];
						$to_print =~s/\/\d$//;
						$to_print =~s/\s\d.+$//;
						delete($reverse_name_list{$to_print});
					}
					
					for(my $i=scalar(@forward_read); $i>0;$i--)
					{
						$deleted_forw++;
						my $to_print = shift(@forward_read);
						print $single_out_handle $to_print;
						my @splitted = split("\n",$to_print);
						my $id = $splitted[0];
						$id =~s/\/\d.*$//;
						$id =~s/\s\d.+$//;
						delete($forward_name_list{$id});
					}
					
					#We flush the matching reads in their respective file
					print $forward_out_handle $read_forw;
					my $to_print = shift(@reverse_read);
					$deleted_rev++;
					print $reverse_out_handle $to_print;
					my @splitted = split("\n",$to_print);
					$to_print = $splitted[0];
					$to_print =~s/\/\d$//;
					$to_print =~s/\s\d.+$//;
					delete($reverse_name_list{$to_print});
					
					#We add the reverse read to the arrays
					$reverse_name_list{$id_rev} = $inserted_rev;
					$inserted_rev++;
					push(@reverse_read, $read_rev);
					
					$found = 1;
				}
				
				if(exists $forward_name_list{$id_rev})
				{
					$pos = $forward_name_list{$id_rev};
					#we flush all the previous reads in the single
					for(my $i=0; $i<$pos-$deleted_forw;$deleted_forw++)
					{
						my $to_print = shift(@forward_read);
						print $single_out_handle $to_print;
						my @splitted = split("\n",$to_print);
						$to_print = $splitted[0];
						$to_print =~s/\/\d.*$//;
						$to_print =~s/\s\d.+$//;
						delete($forward_name_list{$to_print});
					}
					
					#We flush the matching reads in their respective file
					print $reverse_out_handle $read_rev;
					my $to_print = shift(@forward_read);
					$deleted_forw++;
					print $forward_out_handle $to_print;
					my @splitted = split("\n",$to_print);
					$to_print = $splitted[0];
					$to_print =~s/\/\d.*$//;
					$to_print =~s/\s\d.+$//;
					delete($forward_name_list{$to_print});
					
					#We add the forward read to the arrays
					$forward_name_list{$id_forw} = $inserted_forw;
					$inserted_forw++;
					push(@forward_read, $read_forw);
					
					$found = 1;
				}
				
				#If we did not find a match 
				if(!$found)
				{
					$forward_name_list{$id_forw} = $inserted_forw;
					$inserted_forw++;
					push(@forward_read, $read_forw);
					$reverse_name_list{$id_rev} = $inserted_rev;
					$inserted_rev++;
					push(@reverse_read, $read_rev);
					
				}
			}
			
		}
		else
		{
			$end_of_file = 2;
			
			if(exists $reverse_name_list{$id_forw})
			{
				my $pos = $reverse_name_list{$id_forw};
				#we flush all the previous reads in the single
				for(my $i=0; $i<$pos-$deleted_rev;$deleted_rev++)
				{
					my $to_print = shift(@reverse_read);
					print $single_out_handle $to_print;
					my @splitted = split("\n",$to_print);
					$to_print = $splitted[0];
					$to_print =~s/\/\d$//;
					$to_print =~s/\s\d.+$//;
					delete($reverse_name_list{$to_print});
				}
				
				#We flush the matching reads in their respective file
				print $forward_out_handle $read_forw;
				my $to_print = shift(@reverse_read);
				$deleted_rev++;
				print $reverse_out_handle $to_print;
				my @splitted = split("\n",$to_print);
				$to_print = $splitted[0];
				$to_print =~s/\/\d$//;
				$to_print =~s/\s\d.+$//;
				delete($reverse_name_list{$to_print});
			
			}
			else
			{
				print $single_out_handle $read_forw;
			}
		}
	}
	else
	{
		$end_of_file = 1;
	}
}
#If we reached end of forward file
if($end_of_file == 1)
{
	my $id_rev;
	while($id_rev = <$reverse_handle>)
	{
		my $read_rev .= $id_rev;
		chomp($id_rev);
		$id_rev=~s/\/\d$//;
		$id_rev=~s/\s\d.+$//;
		$read_rev .= <$reverse_handle>.<$reverse_handle>.<$reverse_handle>;

		if(exists $forward_name_list{$id_rev})
		{
			my $pos = $forward_name_list{$id_rev};
			#we flush all the previous reads in the single
			for(my $i=0; $i<$pos-$deleted_forw;$deleted_forw++)
			{
				my $to_print = shift(@forward_read);
				print $single_out_handle $to_print;
			}
			
			#We flush the matching reads in their respective file
			print $reverse_out_handle $read_rev;
			my $to_print = shift(@forward_read);
			$deleted_forw++;
			print $forward_out_handle $to_print;
			
		}		
		else
		{
			print $single_out_handle $read_rev;
		}
	}
}
else
{
	my $id_forw;
	while($id_forw = <$forward_handle>)
	{
		my $read_forw .= $id_forw;
		chomp($id_forw);
		$id_forw=~s/\/\d$//;
		$id_forw=~s/\s\d.+$//;
	
		$read_forw .= <$forward_handle>.<$forward_handle>.<$forward_handle>;
		
		if(exists $reverse_name_list{$id_forw})
		{
			my $pos = $reverse_name_list{$id_forw};
			#we flush all the previous reads in the single
			for(my $i=0; $i<$pos-$deleted_rev;$deleted_rev++)
			{
				my $to_print = shift(@reverse_read);
				print $single_out_handle $to_print;
			}
			
			#We flush the matching reads in their respective file
			print $forward_out_handle $read_forw;
			my $to_print = shift(@reverse_read);
			$deleted_rev++;
			print $reverse_out_handle $to_print;
			
		}
		else
		{
			print $single_out_handle $read_forw;
		}
		
	}
}

#We now flush all the remnants reads in the single output
foreach my $to_print (@reverse_read)
{
	print $single_out_handle $to_print;
}

foreach my $to_print (@forward_read)
{
	print $single_out_handle $to_print;
}

my $duration=time() - $start_time;
print("Duration: ". getDuration($duration)." \n");
print("inserted $inserted_forw and deleted $deleted_forw from forward file \n");
print("inserted $inserted_rev and deleted $deleted_rev from reverse file \n");

print (scalar(@forward_read)."\n");
print (scalar(@reverse_read)."\n");
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
        return sprintf('%d minute(s) and %04.2f second(s)', int($duration/60), ($duration%60));
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

	
=pod

=head1 AUTHORS

Gautier SARAH (INRA), gautier.sarah@supagro.inra.fr

=head1 DATE

15.05.2015

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
