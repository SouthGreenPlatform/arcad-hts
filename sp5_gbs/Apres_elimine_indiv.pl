#!/bin/env perl


=pod

=head1 NAME

Apres_elimine_indiv.pl - checks if there are still SNPs that do not segregate and checks
if ,in the descent, there are an individual with genotypic class presents only once, if there 
is the case transforms the genotype into missing data.

=head1 DESCRIPTION

checks if there are still SNPs that do not segregate and checks
if ,in the descent, there are an individual with genotypic class presents only once, if there 
is the case transforms the genotype into missing data

=head1 SYNOPSIS / USAGE

Apres_elimine_indiv.pl	[--help] [--man] [--version]

[-i tabulated-file] \

[-o prefix_output] \

example of input file :

Chromosome	POS	REF	ALT	female	male	Vd011	Vd012	Vd013	Vd015	...

chr1	192497	G	A	AG	GG	GG	GG	GG	AG	...

chr1	235688	G	A	AG	GG	GG	NN	GG	AG	...
...


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
Tabulated file


=item B<[-o]> ([prefix_output]): 
prefix for output file


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

my $out = $prefix_output.".csv";
my $sup = $prefix_output."_removed.csv";

if(! -e $input){ print "the file $input does not exist\n";}

open (F, $input) or die "cannot open the file $input\n";
open(OUT, ">$out") or die"cannot open the file $out\n";
open(SUP, ">$sup") or die"cannot open the file $sup\n";


my (@tab, @tab_base_alt, %CLASSES, $male, $femelle, $classe, $nbr_classes, $chr, $pos, $base_ref, $base_alt, $allele1, $allele2, $alt_base1, $alt_base2, $alt_base3, $b1, $b2, $nophas, $NN, %hash, $taille);

$nophas = 0;
$NN =0;
$nbr_classes =0;
my %count = ();
#my $head = <F>;
while(my $line =<F>){
	chomp $line;
	@tab = split (/\t/, $line);
	$chr = $tab[0];
	$pos = $tab[1];
	$base_ref = $tab[2];
	$base_alt = $tab[3];
	$femelle = $tab[4];
	$male = $tab[5];
	if($line =~ m/^Chromosome/){
		@tab = split (/\t/, $line);
		for(my $i = 0;$i<@tab;$i++){
			print OUT "$tab[$i]\t";
		}
		print OUT "\n";
	}
	else{
	#count genotypic classes
	for(my $i = 6;$i<@tab;$i++){
		$tab[$i] =~ m/(\w)(\w)/;
		# write genotypes in alphabetical order
		my $b1 = $1;
		my $b2 = $2;
		my $GN = "$b1$b2";
		my @GN = split ('', $GN);
		my @ch = sort({$a cmp $b} @GN );
		my $GNT=join('', @ch);
		$classe = "$GNT";
		if((!exists($CLASSES{$classe})) and ($classe ne "NN")){
			$CLASSES{$classe} = 1;
			$nbr_classes++;
		}
		if(exists($CLASSES{$classe}) and ($classe ne "NN")){
			$count{$classe}++;
		}
	}
	# keep SNPs that have more than one genotypic class in the offspring (remove the case: one genotypic class)
	if ($nbr_classes > 1){
		print OUT "$chr\t$pos\t$base_ref\t$base_alt\t$femelle\t$male\t";
		for(my $i = 6;$i<@tab;$i++){
			$tab[$i] =~ m/(\w)(\w)/;
			my $b1 = $1;
			my $b2 = $2;
			my $chaine = "$b1$b2";
			# If a genotypic class is present only once, transform it into NN (missing data)
			my $GN = "$1$2";
			my @GN = split ('', $GN);
			my @ch = sort({$a cmp $b} @GN );
			my $GNT=join('', @ch);
			$classe = "$GNT";
			if($count{"$classe"} == 1){
				$chaine = "NN";
				print SUP "UneClassDiff		$line\n";
			}
			print OUT "$chaine\t";
		}
			print OUT "\n";
	}
	else{
		print SUP "UneClasse	$line\n";
	}
	}
%CLASSES = ();
$nbr_classes = 0;
%count = ();

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

=head1 TODO

=cut
