
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

Cons.pm - All constants variables

=head1 DESCRIPTION

All constants variables

=head1 SYNOPSIS / USAGE

use strict;
use Carp;
use Modules::Config::Cons;

my $gatkPath = &$Softwares::GATK_DIRECTORY or confess($!);

=cut

package Cons;

use strict;
use warnings;

=head1 CONSTANTS

C<$PICARD_TOOLS_DIRECTORY> Path where to find all Picard tools (*.jar), Version 1.7

C<$GATK_DIRECTORY> Path where to find all GATK jar files, Version 1.6

C<$JAVA_PATH> Lastest version of java command, Version 1.6

C<$SAMTOOLS_PATH> Lastest version of samtools, Version 0.1.18

C<$BWA_PATH> Lastest version of bwa, Version 0.5.9

=cut

our $QUEUE = 'arcad.q';
our $QSUB = (
	'-b' => y,
	'-cwd' => '',
	'-q' => $QUEUE,
);

sub qsub
{
	my $command = shift;
	my $job_name = shift;
	
	
}

=head1 INCOMPATIBILITIES

Fully compatble with all version of perl

=head1 BUGS AND LIMITATIONS

NO BUG FOUND YET

=head1 AUTHORS

=over 2

=item Homa Felix

=item felix.homa-at-cirad.fr

=back

=head1 VERSION

1.0.0

=head1 DATE

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
