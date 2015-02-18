
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

Files - File manipulation module

=head1 DESCRIPTION

This modules aimes to store all function about files manipulation, but it dont 
inherit of perl File modules.

=head1 SYNOPSIS / USAGE

[% Not yet commented %]

=cut

package Files;

use strict;
use warnings;
use Error qw(:try);

=head1 METHODS

=head2 getFiles

=head3 description

This fonction return a reference to an array containing all files 
matching with the given pattern

=head3 args

B<input> folder in wich to apply search

B<pattern> search pattren, for example : '*.fastq'

B<subfolder> boolean to specify, if search would be apply to subdirectories

=head3 return

B<Array ref>

=cut
my @a_files;
sub getFiles
{
	my $o_self = shift;
	my $input = shift;
	my $pattern = shift;
	my $subfolder = shift;
	
	throw Error::Simple("No folder specified") if( !$input or $input eq '' );
	$input =~ s/\/*$//;
	
	my $list_command;
	if($subfolder){
		$list_command = "find -L $input -type f -name '$pattern' ";
	}
	else{
		$list_command = "ls -1 $input/$pattern ";
	}
	
	$list_command = `printf \"\%s\n\" \$($list_command)`;

	@a_files = (@a_files, split(/\n/, $list_command));
	
	return \@a_files;
}

=head2 compressFiles

=head3 description

This function zip and encrypt an array of file after checking file existence.
Zipped files have the same name than the original file with extension '.zip'.
A password is required to unzip files.
Throw an exception if a system cammand return an exit code different to 0.

=head3 args

B<files> array of string, string contains all path toward files.

B<keepOriginal> boolean, if false only the zip file is keep on disk

=head3 return

B<zipped files> Array ref, containing all zipped file path.

=cut

sub compressFiles
{
	my $o_self = shift;
	my $ra_files = shift;
	my $keepOriginal = shift;
	my @a_zipped = ();
	
	$keepOriginal = 0 unless( $keepOriginal );
	
	foreach my $file (@$ra_files)
	{
		if( -e $file )
		{
			my $zipname = $file . '.zip';
			my $zipcommand = "zip -9 -r -P \@rc\@d $zipname $file";
			system("$zipcommand");
			throw Error::Simple("unable to zip file $file\n$!\nCODE $?")  if( $? != 0 );
			push @a_zipped, $file;
			unlink $file if( !$keepOriginal  && -e $zipname );
		}
	}
	return \@a_zipped;
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
