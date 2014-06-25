#!/bin/env perl


=pod

=head1 NAME

abTOconsensus.pl - Codage for consensus map. Coding compatible with Carthagene tool

=head1 DESCRIPTION

Transforms abcd coding into consensus coding for Carthagène tool

=head1 SYNOPSIS / USAGE

abTOconsensus.pl	[--help] [--man] [--version]
[-i abcd tabulated-file] \
[-o prefix_output] \

example of input file :

Chromosome	POS	Marker	phase	Segregation	Child1	Child2	Child3 

chr1	192497	chr1_192497	{-0}	<abxaa>	aa	aa	aa

chr1	235688	chr1_235688	{-1}	<abxaa>	aa	--	aa

chr1	778202	chr1_778202	{01}	<abxaa>	ab	ab	ab


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
abcd input file with known phase

=item B<[-o]> ([prefix_output]): 
prefix_output 


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
if ($version) {print "1.0.0\n";}

my $output = "$prefix_output".".txt";
#my $out1 = "$prefix_output"."_male.txt";


if(defined($input) && ! -e $input )
{
	warn "$input not exists ";
	exit(0);
	
}

if(defined($output) && -e $$output )
{
	warn "$$output already exists ";
	exit(0);
	
}

use strict;


open (F, $input)or die "cannot open the file $input\n";
open(OUT, ">$output") or die"cannot open the file $output\n";

my (@tab,$codeF, $codeM, $j, $contig, $marker, $nbr_indiv, $nbr_marker, $nbr_markerM, $F, $M, $marker_miror, $pos, $femelle, $male, $classe, $segregation);
$j = 0;

my $head = <F>;

while(my $line =<F>){
	chomp $line;
	my @head = split (/\t/, $head);
	my @tab1 = split (/\t/, $line);
	my $segregation = $tab1[4];
	my $phase = $tab1[3];
	$nbr_indiv = $#head-4;
	$nbr_marker = 0;
	my $code;
	
 # consensus coding
 # F normal, M normal
	if($phase eq "{00}"){
		if($segregation eq "<abxaa>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "c";
				}
				elsif($tab[$i] eq "aa"){
					$code = "3";
				}
			}
		}
		elsif($segregation eq "<aaxab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "aa"){
					$code = "5";
				}
				elsif($tab[$i] eq "ab"){
					$code = "a";
				}
			}
		}
		elsif($segregation eq "<abxa0>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "4";
				}
				elsif($tab[$i] eq "a-"){
					$code = "3";
				}
				elsif($tab[$i] eq "b0"){
					$code = "8";
				}
			}
		}
		elsif($segregation eq "<a0xab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "2";
				}
				elsif($tab[$i] eq "a-"){
					$code = "5";
				}
				elsif($tab[$i] eq "b0"){
					$code = "8";
				}
			}
		}
		elsif($segregation eq "<abxab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "aa"){
					$code = "1";
				}
				elsif($tab[$i] eq "ab"){
					$code = "6";
				}
				elsif($tab[$i] eq "bb"){
					$code = "8";
				}
			}
		}
		elsif($segregation eq "<abxac>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "4";
				}
				elsif($tab[$i] eq "bc"){
					$code = "8";
				}
				elsif($tab[$i] eq "ac"){
					$code = "2";
				}
				elsif($tab[$i] eq "aa"){
					$code = "1";
				}
			}
		}
		elsif($segregation eq "<abxcd>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ac"){
					$code = "1";
				}
				elsif($tab[$i] eq "ad"){
					$code = "2";
				}
				elsif($tab[$i] eq "bc"){
					$code = "4";
				}
				elsif($tab[$i] eq "bd"){
					$code = "8";
				}
			}
		}
	}
 # F normal, M mirror
	elsif($phase eq "{01}"){
		if($segregation eq "<abxaa>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "c";
				}
				elsif($tab[$i] eq "aa"){
					$code = "3";
				}
			}
		}
		elsif($segregation eq "<aaxab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "aa"){
					$code = "a";
				}
				elsif($tab[$i] eq "ab"){
					$code = "5";
				}
			}
		}
		elsif($segregation eq "<abxa0>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "8";
				}
				elsif($tab[$i] eq "a-"){
					$code = "3";
				}
				elsif($tab[$i] eq "b0"){
					$code = "4";
				}
			}
		}
		elsif($segregation eq "<a0xab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "1";
				}
				elsif($tab[$i] eq "a-"){
					$code = "a";
				}
				elsif($tab[$i] eq "b0"){
					$code = "4";
				}
			}
		}
		elsif($segregation eq "<abxab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "aa"){
					$code = "2";
				}
				elsif($tab[$i] eq "ab"){
					$code = "9";
				}
				elsif($tab[$i] eq "bb"){
					$code = "4";
				}
			}
		}
		elsif($segregation eq "<abxac>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "8";
				}
				elsif($tab[$i] eq "bc"){
					$code = "4";
				}
				elsif($tab[$i] eq "ac"){
					$code = "2";
				}
				elsif($tab[$i] eq "aa"){
					$code = "1";
				}
			}
		}
		elsif($segregation eq "<abxcd>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ac"){
					$code = "2";
				}
				elsif($tab[$i] eq "ad"){
					$code = "1";
				}
				elsif($tab[$i] eq "bc"){
					$code = "8";
				}
				elsif($tab[$i] eq "bd"){
					$code = "4";
				}
			}
		}
	}
 # F mirror, M normal
	elsif($phase eq "{10}"){
		if($segregation eq "<abxaa>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "3";
				}
				elsif($tab[$i] eq "aa"){
					$code = "c";
				}
			}
		}
		elsif($segregation eq "<aaxab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "aa"){
					$code = "5";
				}
				elsif($tab[$i] eq "ab"){
					$code = "a";
				}
			}
		}
		elsif($segregation eq "<abxa0>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "1";
				}
				elsif($tab[$i] eq "a-"){
					$code = "c";
				}
				elsif($tab[$i] eq "b0"){
					$code = "2";
				}
			}
		}
		elsif($segregation eq "<a0xab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "8";
				}
				elsif($tab[$i] eq "a-"){
					$code = "5";
				}
				elsif($tab[$i] eq "b0"){
					$code = "2";
				}
			}
		}
		elsif($segregation eq "<abxab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "aa"){
					$code = "4";
				}
				elsif($tab[$i] eq "ab"){
					$code = "9";
				}
				elsif($tab[$i] eq "bb"){
					$code = "2";
				}
			}
		}
		elsif($segregation eq "<abxac>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "1";
				}
				elsif($tab[$i] eq "bc"){
					$code = "2";
				}
				elsif($tab[$i] eq "ac"){
					$code = "8";
				}
				elsif($tab[$i] eq "aa"){
					$code = "4";
				}
			}
		}
		elsif($segregation eq "<abxcd>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ac"){
					$code = "4";
				}
				elsif($tab[$i] eq "ad"){
					$code = "8";
				}
				elsif($tab[$i] eq "bc"){
					$code = "1";
				}
				elsif($tab[$i] eq "bd"){
					$code = "2";
				}
			}
		}
	}
 # F mirror, M mirror
	elsif($phase eq "{11}"){
		if($segregation eq "<abxaa>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "3";
				}
				elsif($tab[$i] eq "aa"){
					$code = "c";
				}
			}
		}
		elsif($segregation eq "<aaxab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "aa"){
					$code = "a";
				}
				elsif($tab[$i] eq "ab"){
					$code = "5";
				}
			}
		}
		elsif($segregation eq "<abxa0>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "2";
				}
				elsif($tab[$i] eq "a-"){
					$code = "c";
				}
				elsif($tab[$i] eq "b0"){
					$code = "1";
				}
			}
		}
		elsif($segregation eq "<a0xab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "4";
				}
				elsif($tab[$i] eq "a-"){
					$code = "1";
				}
				elsif($tab[$i] eq "b0"){
					$code = "a";
				}
			}
		}
		elsif($segregation eq "<abxab>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "aa"){
					$code = "8";
				}
				elsif($tab[$i] eq "ab"){
					$code = "6";
				}
				elsif($tab[$i] eq "bb"){
					$code = "1";
				}
			}
		}
		elsif($segregation eq "<abxac>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ab"){
					$code = "2";
				}
				elsif($tab[$i] eq "bc"){
					$code = "1";
				}
				elsif($tab[$i] eq "ac"){
					$code = "4";
				}
				elsif($tab[$i] eq "aa"){
					$code = "8";
				}
			}
		}
		elsif($segregation eq "<abxcd>"){
			for(my $i = 4;$i<@head-1;$i++){
				if($tab[$i] eq "ac"){
					$code = "8";
				}
				elsif($tab[$i] eq "ad"){
					$code = "4";
				}
				elsif($tab[$i] eq "bc"){
					$code = "2";
				}
				elsif($tab[$i] eq "bd"){
					$code = "1";
				}
			}
		}
	}

print OUT "data type f2 intercross\n";
print OUT "$nbr_indiv $nbr_marker 0 0\n";

}

close F;
close OUT;



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
