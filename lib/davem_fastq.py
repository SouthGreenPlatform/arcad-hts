# -*- coding: utf-8 -*-

"""
AUTEUR
       Écrit par Vincent MAILLOL

SIGNALER DES BOGUES
       maillol@supagro.inra.fr

COPYRIGHT
       Copyright © 2011 DAVEM.  Licence  GPLv3+ :  GNU
       GPL version 3 ou supérieures <http://gnu.org/licenses/gpl.html>
       Ce  programme  est  un  logiciel  libre.  Vous pouvez le modifier et le
       redistribuer. Il n'y a AUCUNE GARANTIE dans la mesure autorisée par  la
       loi.
"""

import sys

__version__ = "0.0.2"

class Fastq_read( object ) :
    def __init__(self, str_reads ):
        """
        str_read - @name\nsequence\nplus\nquality\n
                    final \n is optional.
        """
        ( self.name,
          self.seq,
          self.plus_line,
          self.qual ) = str_reads.splitlines()

    def cut_start( self, size ) :
        """
        supprime le debut de la  lecture
        """
        if size >= len( self.seq ) :
            raise Exception( "can't cut %d bases. The read is too short" % (size) )
        if size > 0 :
            self.seq  =  self.seq[ size : ]
            self.qual = self.qual[ size : ]

    def get_member( self, format ) :
        """
        return the member of a pair
        format - "Illumina" or "Casava1.8"
        """
        if  format == "Casava1.8" :
            member = self.name.split( None, 1 )[1][0]
            if member == "1" :
                return 1
            elif  member == "2" :
                return 2

        elif format == "Illumina" :
            member = self.name[-1]
            if member == "1" :
                return 1
            elif  member == "2" :
                return 2
        
        else :
            raise ValueError( "format must be 'Illumina' or 'Casava1.8'" ) 

    def cut_end( self, size ) :
        """
        supprime la fin de la  lecture
        """
        if size >= len( self.seq ) :
            raise Exception( "can't cut %d bases. The read is too short" % (size) )
        if size > 0 :
            self.seq  =  self.seq[ : - size ]
            self.qual = self.qual[ : - size ]

    def convert_qual( format ) :
        raise NotImplementedError( "for next time" )

    def get_name_meta_data( self ) :
        """
        return meta data in name read.
        """
        name_and_data = self.name.split( None, 1 )
        if len( name_and_data ) > 1 :
            return name_and_data[1]
        return ""
            
    def split_name( self ) :
        """
        retourne le nom de la séquence sans le prefixe '@'
        et sans les meta donnée qui sont retournées.
        
        return - [cleaned_name, meta_data ]
        """
        return self.name[1:].split( None, 1 )

    def __str__( self ) :
        return "\n".join( ( self.name, self.seq, self.plus_line, self.qual ) )

    def __repr__( self ) :
        return "<fastq read '%s' at 0x%x>" % ( self.name, id( self ) )


class Fastq_file( file ) :
    """
    Pour manipuler les fichiers fastq.
    """
    def __init__( self, path, mode ) :
       file.__init__( self, path, mode )
       self.seq_already_write = False

    def __iter__( self ) :
        return self
 
    def readline( self ) :
        """
        extraire la sequence suivante du fichier.
        """
        r =  super( Fastq_file, self ).readline()
        r += super( Fastq_file, self ).readline()
        r += super( Fastq_file, self ).readline()
        r += super( Fastq_file, self ).readline()
        return r


    def next( self ) :
        r =  super( Fastq_file, self ).next()
        r += super( Fastq_file, self ).next()
        r += super( Fastq_file, self ).next()
        r += super( Fastq_file, self ).next()
        return r

    def write( self, seq ) :
        r"""
        seq doit etre au format @ref\nACTG\n+\nffff
        sans aucun autre \n
        """
        if self.tell() == 0 :
           super( Fastq_file, self ).write( seq )
        else :
           super( Fastq_file, self ).write( "\n" + seq )
            
    def sort( self, path ) :
        """
        Crée une copie triée du fichier fastq.
        """
        index = []
         
        seq = self.next()
        while seq :
             name_seq = seq.split( '\n', 1 )[ 0 ]
             index.append( [ name_seq, self.tell() - len( seq ) ] )
             seq = self.next()

        index.sort()
        
        sorted_file = Fastq_file( path, "w" )
        for name, pos in index :
            self.seek( pos )
            sorted_file.write( self.readline().rstrip( '\n' ) )
        sorted_file.close()
        

        
def clean_seq_name( seq_name ) :
    """
    retourne le nom de la séquence sans le prefixe '@'
    et sans les meta donnée qui sont retournées.
    
    return - [cleaned_name, meta_data ]
    """
    return seq_name[1:].split( None, 1 )



    

if __name__ == "__main__" :
    of = open( "sort_test_fatq.fq", "w" )
    of.write(  """@PNst-A#0 (mate 1, length=101)
ATACCCCAGTCCATAAAACCCGAACCCCAAACCCCAATCCCTAAACCCTAAAACCGTAAAGTCCAAAACGCTAACCCCTTAACCCTAAACCCTAAACCCTG
+
GFGGDGGGGFFGGEGGGGGCGGGGFGGEGFGFDGEGEGGGFFFAFGGBGFBFDG=FCDC<EAE@BGDFF5FAEEEBEEBDDEEGFFFEFFFFEF5DE@FCB
@PNst-F#0 (mate 2, length=101)
CCATACCCCAAGCCCTAAAACCCAAACCCCAAACCCTAGTACTCAAACTCCAAGCCCTAAACCCTAAAACCCTAAAATTCAAAAACCCTAAAATCTAAACC
+
HHHHHHHHGHHEHHHHHHHHHHGHHHHHDHHHFHFHGHHEGGFGG.FGDGGI@GDGGE@GBGGGGFHEHHHHHHFHFHFHBB@9DEBEEC1DEEEFFCFBF
@PNst-C#0 (mate 2, length=101)
AAACCCTAAAACCCTAAAATTCAAAAACCCTAAAATCTAAACCCTAAACCGTACACCCTAAACCCTGAACCCAAAACCCTATACCCCAGTCCATAAAACCC
+
FGGG8HHGHHHHHHHHHHHHHHHHHHHHHHGG/GGGGGEIHGHFHEGGGGEFFFFFHHDBGHFFHDGGGEFHGHHDDBGDBEEEEHHHEDEEFB8DACCDH
@PNst-E#0 (mate 1, length=101)
GAAACCCCATACCCCAAGCCCTAAAACCCAAACCCCAAACCCTAGTACTCAAACTCCAAGCCCTAAACCCTAAAACCCTAAAATTCAAAAACCCTAAAATC
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHCHHHHFHFHFHCHHHHHHHEHHHHFFDBBGFGGFCEDEGGBBCFCDCFDEE@DFF
@PNst-B#0 (mate 2, length=101)
AAACCCTGAACCCAAAACCCTATACCCCAGTCCATAAAACCCGAACCCCAAACCCCAATCCCTAAACCCTAAAACCGTAAAGTCCAAAACGCTAACCCCTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHHHHHHFHHHHHHHHHHHHHHHHHHHHGHHHGHHCHHCHEBFEFDFHHHBHHHHFBFHFFEE?
@PNst-D#0 (mate 2, length=100)
AAACCCTAAAACCCTAAAATTCAAAAACCCTAAAATCTAAACCCTAAACCGTACACCCTAAACCCTGAACCCAAAACACTATACCCCAGTCCATAAAACC
+
@<77<66?@1E<EEEF=EDGF?FCE?@2=AEGEGCCD4C@DCD?DEDBEDECEF@FCEDF?FF9D5:,<5A=??8@@'3A@CDA4C@:@7CCCBA<B77=""" )
    of.close()
    
    ff = Fastq_file( "test_fatq.fq", "r" )
    ff.sort( 'test_fatq_sorted.fq' )
    ff.readline() 
    ff.close()
