
#  
#  Copyright 2014 INRA-CIRAD
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

=head1 NAME

Bam - BAM file manipulation module

=head1 DESCRIPTION

This modules aimes to store all functions about BAM file manipulation.
This use GATK and samtools softwares.

=head1 SYNOPSIS / USAGE

[% Not yet commented %]

=cut

package Bam;

use strict;
use warnings;
use Carp;
use Error qw(:try);
use lib '..';
use Modules::Config::Softwares;

#-----------------------------------------------------
# EXECUTABLES path
#-----------------------------------------------------
my $java_opts = 'Xmx5g'; # Pour les très gros fichiers, option inutile pour l'instant
my $JAVA_PATH = &$Softwares::JAVA_PATH or confess("$!");
my $GATK_DIR  = &$Softwares::GATK_DIRECTORY or confess("$!");
my $GATK_COMMAND = "$JAVA_PATH -jar $GATK_DIR/GenomeAnalysisTK.jar";
my $SAMTOOLS_PATH = &$Softwares::SAMTOOLS_PATH or confess("$!");
#-----------------------------------------------------

=head1 METHODS

=head2 Constructors

=cut

sub new 
{
	my ($proto, $ra_bams) = @_;
	my $class = ref($proto) || $proto;
	
	throw Error::Simple "No valid input" if( !defined $ra_bams || 0 >= scalar $ra_bams );
	my $o_self = {bams => $ra_bams};
	
	bless $o_self, $class;
	return $o_self;
}


=head2 index

=head3 description


=head3 args

=head3 return

=cut

sub indexBam
{
	my $o_self = shift;
	
	for( @{$o_self->{bams}} )
	{
		my $bam = $_;
		if( defined $bam && -e $bam )
		{
			my $bai = "${bam}.bai";
			next if( -e $bai );
			if( $ENV{JOB_ID} )
			{
				$? = system( "$SAMTOOLS_PATH index $bam $bai" );
				throw Error::Simple "Bam index creation failed : $!\njob id '$ENV{JOB_ID}'\n" if( !-e $bai && $? != 0 );
			}
			else
			{
				my $rh_qsub = &$Softwares::RH_QSUB;
				$$rh_qsub{N} = 'samtools_index';
				my $qsub = &Softwares::make_qsub_command($rh_qsub);
				$? = system( "$qsub '$SAMTOOLS_PATH index $bam $bai'" );
				throw Error::Simple("Job 'samtools_index' failed, exit code $?") if( !-e $bai && $? != 0 );
			}
		}
	}
}



=head1 DEPENDENCIES

Nothing

=head1 INCOMPATIBILITIES

<Fully compatble with all version of perl :)>

=head1 BUGS AND LIMITATIONS

<NO BUG FOUND YET>

=head1 AUTHORS

=over 1

=item Felix Homa

=item felix.homa-at-cirad.fr

=back

=head1 VERSION

1.0.0

==head1 DATE

26/06/2012

=head1 LICENSE AND COPYRIGHT

  Copyright 2014 INRA-CIRAD
  
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

1;
