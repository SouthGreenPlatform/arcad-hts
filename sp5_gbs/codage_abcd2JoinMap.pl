#!/bin/env perl


=pod

=head1 NAME

codage_abcd2JoinMap.pl - transforms abcd coding into JoinMap coding 

=head1 DESCRIPTION

transforms abcd coding into JoinMap coding 

=head1 SYNOPSIS / USAGE

codage_abcd2JoinMap.pl	[--help] [--man] [--version]

[-i input-file] \

[-o output] \


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

=item B<[-i]> ([input_file]): 
input file with abcd coding

=item B<[-o]> ([output]): 
output file for JoinMap


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

my ($man, $help, $version, $input, $output);
$man = 0;
$help = 0;
$version = 0;

GetOptions(
	"help|?"      => \$help,
	"man"         => \$man,
	"version"     => \$version,
	"input|i=s"   => \$input,
	"output|o=s"  => \$output
) or pod2usage(2);

if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}
if ($version) {print "1.0.0\n"; exit;}

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

open (F, $input) or die "cannot open the file $input\n";
my $out_tmp = $output."_tmp.txt";
open(OUT, ">$out_tmp") or die "cannot open the file $out_tmp\n";

my (@tab,$codeF, $codeM, $j, $Indiv, $contig, $marker, $nbr_indiv, $nbr_markerF, $nbr_markerM, $F, $M, $marker_miror, $pos, $femelle, $male, $classe, $segregation);
$j = 0;

my $head = <F>;
print OUT "Markers	segregation		individual	";
	my @head = split (/\t/, $head);
	for(my $i=4; $i<@head;$i++){
		print OUT "$head[$i]\t";
	}
	print OUT "\n";

while(my $line =<F>){
	chomp $line;
	my @tab = split (/\t/, $line);
	
	my $marker = $tab[2];
	my $pos = $tab[1];
	my $segregation = $tab[3];
	#my $segregation = $tab[5];

	if($segregation eq "<abxaa>"){
		$segregation = "<lmxll>";
		print OUT "$marker\t$segregation\t\t\t";
		for(my $i = 4;$i<@tab;$i++){
			if ($tab[$i] eq "aa"){
				$Indiv = "ll";
			}
			elsif ($tab[$i] eq "ab"){
				$Indiv = "lm";
			}
			elsif ($tab[$i] eq "--"){
				$Indiv = "--";
			}
			print OUT "$Indiv\t";
		}
		print OUT "\n";
	}
	elsif($segregation eq "<aaxab>"){
		$segregation = "<nnxnp>";
		print OUT "$marker\t$segregation\t\t\t";
		for(my $i = 4;$i<@tab;$i++){
			if ($tab[$i] eq "aa"){
				$Indiv = "nn";
			}
			elsif ($tab[$i] eq "ab"){
				$Indiv = "np";
			}
			elsif ($tab[$i] eq "--"){
				$Indiv = "--";
			}
			print OUT "$Indiv\t";
		}
		print OUT "\n";
	}
	elsif($segregation eq "<abxab>"){
		$segregation = "<hkxhk>";
		print OUT "$marker\t$segregation\t\t\t";
		for(my $i = 4;$i<@tab;$i++){
			if ($tab[$i] eq "aa"){
				$Indiv = "hh";
			}
			elsif ($tab[$i] eq "bb"){
				$Indiv = "kk";
			}
			elsif ($tab[$i] eq "ab"){
				$Indiv = "hk";
			}
			elsif ($tab[$i] eq "--"){
				$Indiv = "--";
			}
			print OUT "$Indiv\t";
		}
		print OUT "\n";
	}
	elsif($segregation eq "<abxac>"){
		$segregation = "<efxeg>";
		print OUT "$marker\t$segregation\t\t\t";
		for(my $i = 4;$i<@tab;$i++){
			if ($tab[$i] eq "aa"){
				$Indiv = "ee";
			}
			elsif ($tab[$i] eq "ab"){
				$Indiv = "ef";
			}
			elsif ($tab[$i] eq "ac"){
				$Indiv = "eg";
			}
			elsif ($tab[$i] eq "bc"){
				$Indiv = "fg";
			}
			elsif ($tab[$i] eq "--"){
				$Indiv = "--";
			}
			print OUT "$Indiv\t";
		}
		print OUT "\n";
	}
	elsif($segregation eq "<abxcd>"){
		$segregation = "<abxcd>";
		print OUT "$marker\t$segregation\t\t\t";
		for(my $i = 4;$i<@tab;$i++){
			if ($tab[$i] eq "ac"){
				$Indiv = "ac";
			}
			elsif ($tab[$i] eq "ad"){
				$Indiv = "ad";
			}
			elsif ($tab[$i] eq "bc"){
				$Indiv = "bc";
			}
			elsif ($tab[$i] eq "bd"){
				$Indiv = "bd";
			}
			elsif ($tab[$i] eq "--"){
				$Indiv = "--";
			}
			print OUT "$Indiv\t";
		}
		print OUT "\n";
	}
}

close F;
close OUT;

TransposerFichier("$out_tmp","$output");

system ("rm $out_tmp" );
#============================================
# TransposerFichier
# Transposer un fichier tabule
#Convertir lignes en colonne
#============================================
sub TransposerFichier {
  my ( $FichierTabuleOriginal, $FichierTranspose ) = @_;
 
  my %HashTranspose;
 
  # Lecture du fichier tabule
  open( my $FH, '<', $FichierTabuleOriginal ) 
    or die("Impossible de lire le fichier $FichierTabuleOriginal\n");
 
  while ( my $Line = <$FH> ) {
    chomp $Line;
    my @data = split( /\t/, $Line );
    for ( my $i = 0; $i < scalar(@data); $i++ ) {
      $HashTranspose{$i} .= $data[$i] . "\t";
    }
  }
  close($FH);
 
  # Création du fichier transpose
  open( my $FHTranpose, '>', $FichierTranspose ) 
    or die("Impossible de creer le fichier $FichierTranspose\n");
 
  foreach ( sort { $a <=> $b } keys %HashTranspose ) {
    $HashTranspose{$_} =~ s{\t$}{};
    print {$FHTranpose} "$HashTranspose{$_}\n";
  }
  close($FHTranpose);
 
  return $FichierTranspose;
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

=head1 TODO

=cut

