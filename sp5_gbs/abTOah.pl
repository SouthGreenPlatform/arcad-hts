#!/bin/env perl


=pod

=head1 NAME

abTOah.pl - Transformes ab coding into AH coding. Coding compatible with Carthagene tool

=head1 DESCRIPTION

Transformes ab coding into AH coding. Coding compatible with Carthagene tool

Creates two files : female genotype file and male genotype file

=head1 SYNOPSIS / USAGE

abTOah.pl	[--help] [--man] [--version]
[-i tabulated-file] \
[-o prefix_output] \


example of input file :

Chromosome	POS	Marker	Segregation	Child1	Child2	Child3 

chr1	192497	chr1_192497	<abxaa>	aa	aa	aa

chr1	235688	chr1_235688	<abxaa>	aa	--	aa

chr1	778202	chr1_778202	<abxaa>	ab	ab	ab

=cut

use strict;
use Pod::Usage;
use Getopt::Long;
use Carp qw (cluck confess croak);


=pod

=head1 OPTIONS

=over 8

=item B<--help>
Show the brief help information.

=item B<--man>
Read the manual, with examples.

=item B<--version>
Show the version number and exit.

=head2 required parameters

=over 5

=item B<[-i]> ([input_file]): 
input tabulated file


=item B<[-o]> ([prefix_output]): 
prefix_output for male and femele output


=cut

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

my ($man, $help, $version, $input, $prefix_output);
$man = 0;
$help = 0;
$version = 0;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"input|i=s"   => \$input,
	"output|o=s"  => \$prefix_output	
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
if ($version) {print "1.0.0\n"; exit}

my $out = "$prefix_output"."_female.txt";
my $out1 = "$prefix_output"."_male.txt";


if(defined($input) && ! -e $input )
{
	warn "$input not exists ";
	exit(0);
	
}

if(defined($out) && -e $out )
{
	warn "$out already exists ";
	exit(0);
	
}

if(defined($out1) && -e $out1 )
{
	warn "$out1 already exists ";
	exit(0);
	
}

open (F, $input)or die "cannot open the file $input\n";
open(OUT, ">$out") or die"cannot open the file $out\n";
open(OUT1, ">$out1") or die"cannot open the file $out1\n";

my (@tab,$codeF, $codeM, $j, $contig, $marker, $nbr_indiv, $nbr_markerF, $nbr_markerM, $F, $M, $marker_miror, $pos, $femelle, $male, $classe, $segregation);
$j = 0;

my $head = <F>;
while(my $line =<F>){
	chomp $line;
	my @head = split (/\t/, $head);
	my @tab1 = split (/\t/, $line);
	my $segregation = $tab1[3];
	$nbr_indiv = $#head-4;
	
	if($segregation eq "<abxaa>"){
	$F++;
	}
	elsif($segregation eq "<aaxab>"){
		$M++;
	}
	elsif($segregation eq "<a0xab>"){
		$F++;
		$M++;
	}
	elsif($segregation eq "<abxa0>"){
		$F++;
		$M++;
	}
	elsif($segregation eq "<abxab>"){
		$F++;
		$M++;
	}
	elsif($segregation eq "<abxac>"){
		$F++;
		$M++;
	}
	elsif($segregation eq "<abxcd>"){
		$F++;
		$M++;
	}
}
$nbr_markerF = $F*2;
$nbr_markerM = $M*2;

print OUT "data type f2 backcross\n";
print OUT "$nbr_indiv $nbr_markerF 0 0\n";

print OUT1 "data type f2 backcross\n";
print OUT1 "$nbr_indiv $nbr_markerM 0 0\n";

close F;
open (F, $input)or die "cannot open the file $input\n";

while(my $line =<F>){
	chomp $line;
	@tab = split (/\t/, $line);
	my @head = split (/\t/, $head);

	$contig = $tab[0];
	$pos = $tab[1];
	$marker = $tab[2];
	$marker_miror = "$marker"."_m";
	$segregation = $tab[3];
	
	if($segregation eq "<abxaa>"){
		print OUT "*$marker ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "aa"){
				$codeF = "A";
			}
			elsif($classe eq "ab"){
				$codeF = "H";
			}
			elsif($classe eq "--"){
				$codeF = "-";
			}
			print OUT "$codeF ";
		}
		print OUT "\n";
		print OUT "*$marker_miror ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "aa"){
				$codeF = "H";
			}
			elsif($classe eq "ab"){
				$codeF = "A";
			}
			elsif($classe eq "--"){
				$codeF = "-";
			}
			print OUT "$codeF ";
		}
		print OUT "\n";
	}
	elsif($segregation eq "<aaxab>"){
		print OUT1 "*$marker ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "aa"){
				$codeM = "A";
			}
			elsif($classe eq "ab"){
				$codeM = "H";
			}
			elsif($classe eq "--"){
				$codeM = "-";
			}
			print OUT1 "$codeM ";
		}
		print OUT1 "\n";
		print OUT1 "*$marker_miror ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "aa"){
				$codeM = "H";
			}
			elsif($classe eq "ab"){
				$codeM = "A";
			}
			elsif($classe eq "--"){
				$codeM = "-";
			}
			print OUT1 "$codeM ";
		}
		print OUT1 "\n";
	}
	elsif($segregation eq "<a0xab>"){
		print OUT "*$marker ";
		print OUT1 "*$marker ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "a-"){
				$codeF = "-";
				$codeM = "A";
			}
			elsif($classe eq "ab"){
				$codeF = "A";
				$codeM = "H";
			}
			elsif($classe eq "b0"){
				$codeF = "H";
				$codeM = "H";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
		print OUT "*$marker_miror ";
		print OUT1 "*$marker_miror ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "a-"){
				$codeF = "-";
				$codeM = "H";
			}
			elsif($classe eq "ab"){
				$codeF = "H";
				$codeM = "A";
			}
			elsif($classe eq "b0"){
				$codeF = "A";
				$codeM = "A";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
	}
	elsif($segregation eq "<abxa0>"){
		print OUT "*$marker ";
		print OUT1 "*$marker ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "a-"){
				$codeF = "A";
				$codeM = "-";
			}
			elsif($classe eq "ab"){
				$codeF = "H";
				$codeM = "A";
			}
			elsif($classe eq "b0"){
				$codeF = "H";
				$codeM = "H";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
		print OUT "*$marker_miror ";
		print OUT1 "*$marker_miror ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "a-"){
				$codeF = "H";
				$codeM = "-";
			}
			elsif($classe eq "ab"){
				$codeF = "A";
				$codeM = "H";
			}
			elsif($classe eq "b0"){
				$codeF = "A";
				$codeM = "A";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
	}
	elsif($segregation eq "<abxab>"){
		print OUT "*$marker ";
		print OUT1 "*$marker ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "aa"){
				$codeF = "A";
				$codeM = "A";
			}
			elsif($classe eq "ab"){
				$codeF = "-";
				$codeM = "-";
			}
			elsif($classe eq "bb"){
				$codeF = "H";
				$codeM = "H";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
		print OUT "*$marker_miror ";
		print OUT1 "*$marker_miror ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "aa"){
				$codeF = "H";
				$codeM = "H";
			}
			elsif($classe eq "ab"){
				$codeF = "-";
				$codeM = "-";
			}
			elsif($classe eq "bb"){
				$codeF = "A";
				$codeM = "A";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
	}
	elsif($segregation eq "<abxac>"){
		print OUT "*$marker ";
		print OUT1 "*$marker ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "aa"){
				$codeF = "A";
				$codeM = "A";
			}
			elsif($classe eq "ab"){
				$codeF = "H";
				$codeM = "A";
			}
			elsif($classe eq "bc"){
				$codeF = "H";
				$codeM = "H";
			}
			elsif($classe eq "ac"){
				$codeF = "A";
				$codeM = "H";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
		print OUT "*$marker_miror ";
		print OUT1 "*$marker_miror ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "aa"){
				$codeF = "H";
				$codeM = "H";
			}
			elsif($classe eq "ab"){
				$codeF = "A";
				$codeM = "H";
			}
			elsif($classe eq "bc"){
				$codeF = "A";
				$codeM = "A";
			}
			elsif($classe eq "ac"){
				$codeF = "H";
				$codeM = "A";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
	}
	elsif($segregation eq "<abxcd>"){
		print OUT "*$marker ";
		print OUT1 "*$marker ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "ac"){
				$codeF = "A";
				$codeM = "A";
			}
			elsif($classe eq "ad"){
				$codeF = "A";
				$codeM = "H";
			}
			elsif($classe eq "bc"){
				$codeF = "H";
				$codeM = "A";
			}
			elsif($classe eq "bd"){
				$codeF = "H";
				$codeM = "H";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
		print OUT "*$marker_miror ";
		print OUT1 "*$marker_miror ";
		for(my $i = 4;$i<@head-1;$i++){
			$classe = $tab[$i];
			if($classe eq "ac"){
				$codeF = "H";
				$codeM = "H";
			}
			elsif($classe eq "ad"){
				$codeF = "H";
				$codeM = "A";
			}
			elsif($classe eq "bc"){
				$codeF = "A";
				$codeM = "H";
			}
			elsif($classe eq "bd"){
				$codeF = "A";
				$codeM = "A";
			}
			elsif($classe eq "--"){
				$codeF = "-";
				$codeM = "-";
			}
			print OUT "$codeF ";
			print OUT1 "$codeM ";
		}
		print OUT "\n";
		print OUT1 "\n";
	}
}

close F;
close OUT;
close OUT1;

=head1 INCOMPATIBILITIES

Fully compatble with any perl version

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 2

=item Hajar CHOUIKI (INRA), hajar.chouiki@supagro.inra.fr

=back

=head1 VERSION

1.0.0

=head1 DATE

[% DATE %]

=head1 LICENSE AND COPYRIGHT

<Please cite Cirad AGAP Data Integration team, if you use this module.
Thank you!> 

=cut
