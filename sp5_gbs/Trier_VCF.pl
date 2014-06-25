#!/bin/env perl


=pod

=head1 NAME

Trier_VCF.pl - Generates a vcf file with sorted SNPs, SNPs that are in the text file

=head1 DESCRIPTION

Generates a vcf file with sorted SNPs, SNPs that are in the text file

=head1 SYNOPSIS / USAGE

Trier_VCF.pl	[--help] [--man] [--version]
[-vcf vcf-file] \
[-i input-file] \
[-o output] \

example of input file (-i) :

Chromosome	POS	Marker	REF	ALT	Segregation	Child1	Child2	Child3

chr1	192497	chr1_192497	<abxaa>	aa	aa	aa

chr1	235611	chr1_235611	<abxaa>	aa	--	aa

chr1	235683	chr1_235683	<abxaa>	aa	--	aa

=cut

use strict;
use Pod::Usage;
use Getopt::Long;
use Carp qw (cluck confess croak);
use Statistics::Distributions ; 

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

=item B<[-vcf]> ([vcf_file]): 
VCF file (VCF file before filters)

=item B<[-i]> ([input_file]): 
input file (output file "eliminer_distortion_kideux.pl" script )

=item B<[-o]> ([output]): 
output file



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

my ($man, $help, $version, $input, $output, $vcf);
$man = 0;
$help = 0;
$version = 0;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"vcf|vcf=s"   => \$vcf,
	"input|i=s"   => \$input,
	"output|o=s"  => \$output
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
if ($version) {print "1.0.0\n"; exit;}

if($vcf eq '')
{
	warn "you forgot the -vcf parameter ";
	exit;
	
}
if(defined($vcf) && ! -e $vcf )
{
	warn "$vcf not exists ";
	exit(0);
	
}

if(defined($input) && ! -e $input )
{
	warn "$input not exists ";
	exit(0);
	
}
if(not defined($output))
{
	warn "you forgot the -o parameter ";
	exit(0);
	
}
if(defined($output) && -e $output )
{
	warn "$output already exists ";
	exit(0);
	
}

open (F1, $vcf) or die "cannot open the file $vcf\n";
open (F2, $input) or die "cannot open the file $input\n";
open(OUT, ">$output") or die"cannot open the file $output\n";

my (@tab, @POS);

my $head = <F2>;
while(my $line =<F2>){
	chomp $line;
	@tab = split (/\t/, $line);
	my $marker = $tab[2];
	$marker =~s /_i//; 
	push(@POS, $marker);
}
close F2;

my @POS1=doublons_tranche(\@POS);

while(my $line =<F1>){
	chomp $line;
	if ($line =~ m/^#/) {
		print OUT "$line\n";
	}
	else{
		my @col = split /\t/, $line;
		my $mark =$col[0]."_".$col[1];
		foreach my $val(@POS1){
			if ($mark eq $val){
				print OUT "$line\n";
			}
		}
	}
}

close F1;
close OUT;

sub doublons_tranche {
	my ($ref_tabeau) = @_;
	my %hash_sans_doublon;
	@hash_sans_doublon{ @{$ref_tabeau} } = ();
	return keys %hash_sans_doublon;
}


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
