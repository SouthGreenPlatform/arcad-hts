
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

ArcadSge.pm - Qsub manager

=head1 DESCRIPTION

[% to fill %]

=head1 SYNOPSIS / USAGE

[% Not yet commented %]

=cut

package ArcadSge;

use strict;
use warnings;
use Error qw(:try);
use Schedule::SGE;

our @ISA = ("Schedule::SGE");

=head1 METHODS

=head2 Constructors

=cut

sub new 
{
	my ($proto, $rh_params) = @_;
	my $class = ref($proto) || $proto;
	
	my $project = $$rh_params{project} if( defined $$rh_params{project} );
	my $mailto  = $$rh_params{mailto}  if( defined $$rh_params{mailto}  );
	my $verbose = $$rh_params{verbose} if( defined $$rh_params{verbose} );
	my $o_self = $class->SUPER::new(
		-project    => $project,
		-mailto     => $mailto,
		-executable =>
		{
			qsub  => '/SGE/n1ge62/bin/lx24-amd64/qsub',
			qstat => '/SGE/n1ge62/bin/lx24-amd64/qstat'
		},
		-verbose    => $verbose
	);
	
	bless $o_self, $class;
	return $o_self;
}

=head2 getFiles

=head3 description

=head3 args

=head3 return

=cut

sub printStatus
{
	my $o_self = shift;
	
	$o_self->user("sarah1");
	my $user = $o_self->user();;
	my %h_status = %{$o_self->status};
	while( my($node, $ra_node_infos) = each( %h_status ) )
	{
		print "Node : $user $node\n";
		if( $node =~ m/bioinfo/ ){
		
		foreach( my($queueType, $nb_proc, $loadAverage, $state) = @$ra_node_infos )
		{
			print "
0. Queue type      = $queueType
1. Processors used = $nb_proc
2. Load average    = $loadAverage
3. State           = $state
";
		}
		}
	}
	
}

sub print_all_jobs
{
	my $o_self = shift;
	
	my @a_status = $o_self->all_jobs;
	foreach my $line ( @a_status )
	{
		print "$line\n";
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

1;
