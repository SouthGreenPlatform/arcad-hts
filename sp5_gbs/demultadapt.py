#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
AUTHOR
       Written by Vincent MAILLOL (modified by Gautier Sarah)

BUGS
       sarah@supagro.inra.fr

COPYRIGHT
       Copyright © 2011 DAVEM, 2014 AGAP.  Licence  GPLv3+ :  GNU
       GPL version 3 ou supérieures <http://gnu.org/licenses/gpl.html>
       This program is free software; you can redistribute it and/or modify
	   it under the terms of the GNU General Public License as published by
	   the Free Software Foundation; either version 3 of the License, or
       (at your option) any later version.
  
       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.
  
       You should have received a copy of the GNU General Public License
       along with this program; if not, see <http://www.gnu.org/licenses/> or 
       write to the Free Software Foundation, Inc., 
       51 Franklin Street, Fifth Floor, Boston,
       MA 02110-1301, USA.
  
  
"""

import sys, os
print os.getcwd()
sys.path.append("")
from davem_fastq import Fastq_read, Fastq_file
import argparse
from itertools import izip
from bisect import bisect_left
 
if sys.version_info[0] == 2:
	if sys.version_info[1] < 6:
		msg = "ERROR: Python should be in version 2.6 or higher"
		sys.stderr.write("%s\n\n" % msg)
		sys.exit(1)
 
class FastqFileType( object ) :
    """
    Fabrique Fastq_file
    """    
    def __init__( self, mode ) :
        self.mode = mode
    
    def __call__( self, path_name) :
        return Fastq_file( path_name, self.mode )


class Selector( object ) :
    """
    Abstract class to look for a line in a table_adaptator.

    table_adaptator is like :
       [ (adaptator-1, output-file-A1, output-file-B1),
         (adaptator-2, output-file-A2, output-file-B2),
                ...
         (adaptator-N, output-file-AN, output-file-BN)
        ]
	
	In single end mode, table adaptator only have one output file by tuple

    You must implement methods __single_select and __paired_select.
    """
    def __init__(self, table_adaptator, single_end) :
        """
        If single_end is True, a call to monSelector.select( sequence )
        will execute the method _single_select_ otherwise, the call will be
        monSelector.select( sequence-1, sequence-2 ) and will execute the 
        method _paired_select_
        
        """
        self.table_adaptator = table_adaptator
        if single_end :
            self.select = self._single_select
        else :
            self.select = self._paired_select

    def _single_select( self, sequence ) :
        """
        Look for a line in table_adaptator with only one sequence
        """
        raise NotImplementedError

    def _paired_select( self, sequence_1, sequence_2 ) :
        """
        Look for a line in table_adaptator with two sequences
        """
        raise NotImplementedError


class Levenshtein_selector( Selector ) :
    table_adaptator = None
    single_end = False
    rate = 0
    def __init__( self, table_adaptator, single_end, rate ) :
        if not isinstance( rate, float ) :
            raise ValueError( "rate argument must be a float not %s" % type( rate ) )
        Selector.__init__( self, table_adaptator, single_end)
        self.rate = rate

    def _single_select( self, sequence) :
        from Levenshtein import ratio
        
        distances = []
        for (adaptator, output_file) in self.table_adaptator :
            dist = ratio( adaptator, sequence[ : len( adaptator ) ] )
            if dist == 1.0 :
                return (adaptator, output_file)
        
            distances.append( dist )
        
        max_dist = max( distances )
        if max_dist >= self.rate and distances.count( max_dist ) == 1 :
            return self.table_adaptator[ distances.index( max_dist ) ]
        
        return None

    def _paired_select( self, sequence_1, sequence_2) :
        from Levenshtein import ratio
        distances_1 = []
        distances_2 = []

        for line in self.table_adaptator :
            adaptator = line[ 0 ]
            dist_1 = ratio( adaptator, sequence_1[ : len( adaptator ) ] )
            dist_2 = ratio( adaptator, sequence_2[ : len( adaptator ) ] )        
            distances_1.append( dist_1 )
            distances_2.append( dist_2 )

        max_dist_1 = max( distances_1 )
        max_dist_2 = max( distances_2 )

        if max_dist_1 > max_dist_2 :
            if max_dist_1 >= self.rate and distances_1.count( max_dist_1 ) == 1 :
                return self.table_adaptator[ distances_1.index( max_dist_1 ) ]
                                
        elif max_dist_1 < max_dist_2 :
            if max_dist_2 >= self.rate and distances_2.count( max_dist_2 ) == 1 :
                return self.table_adaptator[ distances_2.index( max_dist_2 ) ]

        else :
            if max_dist_1 >= self.rate :
                if distances_1.count( max_dist_1 ) == 1 :
                    index_1 = distances_1.index( max_dist_1 )
                    index_2 = distances_2.index( max_dist_2 )
                    if index_1 == index_2 :
                        return self.table_adaptator[ index_1 ]

                elif distances_2.count( max_dist_2 ) == 1 :
                    index_1 = distances_1.index( max_dist_1 )
                    index_2 = distances_2.index( max_dist_2 )
                    if index_1 == index_2 :
                        return self.table_adaptator[ distances_2.index( max_dist_2 ) ]

        return None


class LevenshteinAllSelector( Levenshtein_selector ) :
    """
    Same as Levenshtein_selector except that in paired-end, both members
    of the pair must be above or equal to the min ratio and adaptators of
    both members must be identical
    """

    def _paired_select( self, sequence_1, sequence_2) :
        from Levenshtein import ratio
        distances_1 = []
        distances_2 = []

        for line in self.table_adaptator :
            adaptator = line[ 0 ]
            dist_1 = ratio( adaptator, sequence_1[ : len( adaptator ) ] )
            dist_2 = ratio( adaptator, sequence_2[ : len( adaptator ) ] )        
            distances_1.append( dist_1 )
            distances_2.append( dist_2 )

        max_dist_1 = max( distances_1 )
        max_dist_2 = max( distances_2 )

        if ( max_dist_1 >= self.rate and max_dist_2 >= self.rate 
           and distances_1.count( max_dist_1 ) == distances_2.count( max_dist_2 ) == 1 ) :
               adapt_1 = self.table_adaptator[ distances_1.index( max_dist_1 ) ] 
               adapt_2 = self.table_adaptator[ distances_2.index( max_dist_2 ) ]
               if adapt_1 == adapt_2 :
                    return adapt_1
        else :
            return None

class Std_selector( Selector ):
    """
    Dichotomic search in list_adaptator
    table_adaptator
    If provided index is empty, return None
    """

    def _paired_select( self, sequence_1, sequence_2):
        l1 = self._single_select( sequence_1 )
        l2 = self._single_select( sequence_2 )
        if l1 is None :
            return l2
        
        if l2 is None :
            return l1
        
        if l1 == l2 :
            return l1

        return None
        
        
    def _single_select( self, sequence):
        a = 0
        b = len( self.table_adaptator ) -1
        if b == -1 :
            return None
        
        while a <= b  :
            m = ( a + b ) // 2
            adaptator = self.table_adaptator[ m ][ 0 ]
            start_seq = sequence[ : len( adaptator ) ]
            
            if adaptator > start_seq :
                b = m - 1
            elif adaptator < start_seq :
                a = m + 1
            else :
                return self.table_adaptator[ m ]

        if adaptator == sequence[ : len( adaptator ) ] :
            return self.table_adaptator[ m ]
        return None


def get_adapt_counter( opened_adapt_file ) :
    """
    Return a hash where keys are the adaptators
    and values are initialized with [ name_tag, 0 ]
    """
    d = {}
    opened_adapt_file.seek(0)
    for line in opened_adapt_file :
        if not line.isspace() :
            try :    
                adapt, name_tag = line.split()
            except ValueError :
                print >> sys.stderr, "File '%s' is malformed." %  opened_adapt_file.name
                exit( 1 )
            d[ adapt ] = [ name_tag, 0 ]
    return d


def get_maximal_annalogie( file_adapt ) :
    """
    Compute maximal levenshtein between all adaptators
    """
    from Levenshtein import ratio
    adaptators = []
    for line in file_adapt :
        if line :
            (adapt, name) = line.split()
            if adapt != "*" :           
                adaptators.append( adapt )
   
    ratio_max = 0.0
    for i, adapt in enumerate( adaptators ) :
        for adapt2 in adaptators[i+1:] :
            ratio_max = max( ratio_max,ratio( adapt, adapt2 ) )

    return ratio_max




def get_output_files( opened_adapt_file, prefix, paired_end=True ) :
    """
    Create output files and put them in a list:
	
	if paired_end is True, twa files by adaptator are created
		[ (adaptator, output_file_1, output_file_2 ), ...  ]
		
	otherwise only one
		 [ (adaptator, output_file ), ...  ]
    
    The function return the files table and the Trash files
    Two trash files for paired-end, one for single-end
        
        ( table, (trash-file, ) )
        
    """ 
    

    ada_files = []
    default = None
    cache_name_file_by_adapt = {}

    for line in opened_adapt_file :
        if not line.isspace() :
                try :    
                    adapt, suffix_file = line.split()
                except ValueError :
                    print >> sys.stderr, "File '%s' is malformed." %  opened_adapt_file.name
                    exit( 1 )

                if paired_end :
                    if line[0] == '*' :
                        default = ( Fastq_file( "%s-%s_1.fastq" % (prefix, suffix_file), "w" ),
                                    Fastq_file( "%s-%s_2.fastq" % (prefix, suffix_file), "w" ), )

                    else :
                        if suffix_file in cache_name_file_by_adapt :
                            f1, f2 = cache_name_file_by_adapt[ suffix_file ]
                            ada_files.append( ( adapt, f1, f2 ) )

                        else :
                            f1 = Fastq_file( "%s-%s_1.fastq" % (prefix, suffix_file), "w" )
                            f2 = Fastq_file( "%s-%s_2.fastq" % (prefix, suffix_file), "w" )
                            ada_files.append( (adapt, f1, f2) )
                            cache_name_file_by_adapt[ suffix_file ] = (f1, f2)


                else :
                    # TODO Make cache system for single mode.
                    if line[0] == '*' :
                        default = ( Fastq_file( "%s-%s.fastq" % (prefix, suffix_file), "w" ) , )

                    else :
                        if suffix_file in cache_name_file_by_adapt :
                            f1 = cache_name_file_by_adapt[ suffix_file ]
                            ada_files.append(  ( adapt, f1 ) )
                        else:
                            f1 = Fastq_file( "%s-%s.fastq" % (prefix, suffix_file), "w" )
                            ada_files.append( ( adapt, f1 ) )
                            cache_name_file_by_adapt[ suffix_file ] = ( f1 )
                            
    if default is None :
        print >> sys.stderr, "File '%s' doesn't have a line with the joker tag *.\nAdd a line '*    tag_name_for_trash'." %  opened_adapt_file.name
        sys.exit(1)
    
    ada_files.sort()
    return ada_files, default


def parse_user_argument() :
    """
    Recover user argument
    """
    parser = argparse.ArgumentParser( description="demultiplex fastq_file(s)" )
    
    parser.add_argument( '-V', '--version', action='version', help="Print the version and license",
                         version="%(prog)s 1.1\nCopyright (C) 2011 DAVEM, 2014 AGAP\nGPL3+\nWritten by Vincent Maillol" )

    parser.add_argument( '-v', '--verbose', dest="verbose", action='store_true',
                            help="Be verbose" )
    
    parser.add_argument( 'file_adapt', metavar="FILE_ADAPT", nargs=1, type=argparse.FileType('r'),
                         help="Format is one line by adaptor, such as: adaptor_1<tab>id_sample_1, etc. Last line should be like *<tab>name_trash")

    parser.add_argument( '-f', '--fastq_1', dest="fastq_1", type=FastqFileType( "r" ), action='store', 
                            help="For a single-end file or the first paired-end file" )

    parser.add_argument( '-F', '--fastq_2', dest="fastq_2", type=FastqFileType( "r" ), action='store', default=None,
                            help="For the 2nd paired-end file" )

    parser.add_argument( '-p', '--output_prefix', dest="output_prefix", default="", action='store',
                            help="Output files have name PREFIX-ADAPTOR.fastq"  )

    parser.add_argument( '-l', '--levenshtein', dest="levenshtein", action='store', type=float, default=None,
                            help="Use a Levenshtein distance to demultiplex" )

    parser.add_argument( '-a', '--analogy', dest="analogy", action='store_true',
                            help="Compute the maximal Levenshtein ratio between adaptors" )

    parser.add_argument( '--all', dest="all", action='store_true',
                         help="If this option is used with option levenshtein in paired-end, both members should be higher than the ratio and each should be close to one adaptor. If option levenshtein is not used, this option is not used either." )

    user_args = parser.parse_args()
    user_args.file_adapt = user_args.file_adapt[0]
    user_args.single_end = user_args.fastq_2 is None 
    return user_args

def main() :
    user_args = parse_user_argument()

    if user_args.analogy :
        print "Maximal Levenshtein ratio between adaptors is %f" % get_maximal_annalogie( user_args.file_adapt )
        sys.exit(0)        
            
    output_files_by_adapt, defaults_files = get_output_files( user_args.file_adapt,
                                                              user_args.output_prefix,
                                                              not user_args.single_end )

    nb_reads_writen = get_adapt_counter( user_args.file_adapt )

    user_args.file_adapt.close()

    if user_args.levenshtein : 
        if user_args.all :
            select_output_file = LevenshteinAllSelector( output_files_by_adapt, 
                                                       user_args.single_end,
                                                       user_args.levenshtein )

        else :
            select_output_file = Levenshtein_selector( output_files_by_adapt, 
                                                       user_args.single_end,
                                                       user_args.levenshtein )
    else :
        select_output_file = Std_selector( output_files_by_adapt, 
                                           user_args.single_end )

    if user_args.single_end :
        print "single end"
        default_file = defaults_files[0]
        for str_read in user_args.fastq_1 :
            read = Fastq_read( str_read )
            adapt_and_line = select_output_file.select( read.seq )
            if adapt_and_line is None :
                if user_args.verbose :
                    print "Read '%s' start with %s... and go to *" % (read.name, read.seq[ : 14 ])
                default_file.write( "%s" % str( read ) )
                nb_reads_writen[ '*' ][ 1 ] += 1

            else :
                (adapt, output_file) = adapt_and_line
                if user_args.verbose :
                    print "Read '%s' start with %s... and go to %s" % (read.name, read.seq[ : len( adapt ) ], adapt)
                read.cut_start( len( adapt ) )
                output_file.write( "%s" % str( read ) )
                nb_reads_writen[ adapt ][ 1 ] += 1

        user_args.fastq_1.close()

        for adapt, output_file in output_files_by_adapt :
            if not output_file.closed:
                output_file.write("")
            output_file.close()
        if not default_file.closed:
            default_file.write("")
        default_file.close()

    else :
        print "paired-end"
        (default_file_1, default_file_2) = defaults_files

        for str_read_1, str_read_2 in izip( user_args.fastq_1, user_args.fastq_2 ) :
            read_1 = Fastq_read( str_read_1 )
            read_2 = Fastq_read( str_read_2 )

            adapt_and_line = select_output_file.select( read_1.seq, read_2.seq )
            
            if adapt_and_line is None :
                default_file_1.write( "%s" % str( read_1 ) )
                default_file_2.write( "%s" % str( read_2 ) )
                nb_reads_writen[ '*' ][1] += 1

            else :
                (adapt, output_file_1, output_file_2 ) = adapt_and_line

                read_1.cut_start( len( adapt ) )
                read_2.cut_start( len( adapt ) )

                output_file_1.write( "%s" % str( read_1 ) )
                output_file_2.write( "%s" % str( read_2 ) )
                nb_reads_writen[ adapt ][1] += 1

        user_args.fastq_1.close()
        user_args.fastq_2.close()

        for adapt, file_1, file_2 in output_files_by_adapt :
            if not file_1.closed:
                file_1.write("")
            file_1.close()
            if not file_2.closed:
                file_2.write("")
            file_2.close()
        if not default_file_1.closed:
            default_file_1.write("")
        default_file_1.close()
        if not default_file_2.closed:
            default_file_2.write("")
        default_file_2.close()

    # show stat.
    for nb_reads_by_name in nb_reads_writen.values() :
        print "%s %d reads" % tuple( nb_reads_by_name )


if __name__ == '__main__':
    main()
