from ensembl.seqregion import *
from ensembl.datamodel import *
import MySQLdb


def getDriver(host, user, dbname):
    '''Provides a method to instantiate a Driver object'''
    if Driver._instance == None:
        Driver._instance = Driver(host, user, dbname)
    return Driver._instance


class Driver (object):
    '''Provide access to a database in the core ensembl database system'''
    
    _instance = None

    def __init__(self, host, user, dbname):
        self.conn = MySQLdb.connect(host, user)
        self.db = dbname
        self.tb_adaptor = {'exon': ExonAdaptor, 'gene': GeneAdaptor, 'transcript': TranscriptAdaptor, 'seq_region': SeqregionAdaptor, 'translation': TranslationAdaptor, 'gene_stable_id': GeneStableIdAdaptor, 'transcript_stable_id': TranscriptStableIdAdaptor, 'translation_stable_id': TranslationStableIdAdaptor, 'exon_stable_id': ExonStableIdAdaptor}

    def getAdaptor(self, tbname):    
        adaptor_name = self.tb_adaptor[tbname]
        return adaptor_name(self.db, self.conn.cursor())

    def fetch_sequence_by_region(self, chromosome, start, end, strand):
        '''Obtain the DNA sequence of a particular genomic region (defined by
chromosome, start, end, strand).
        Note: the start and end are based on ensembl 1-offset coordinate system.
        '''
        
        # create a SeqregionAdaptor object
        seq_region = self.getAdaptor('seq_region')
        # obtain a SeqregionAdaptor table object
        seq_regiontb = seq_region.tbobj
        
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        # create a SeqRegion object
        srdb = SeqRegion(seq_regiontb, {17:hg18}, {17:'chr'})
        
        # retrieve the required seq_region_id from the seq_region table by chromosome name
        cursor = seq_region.cursor
        #n = srdb.seqRegionDB.cursor.execute('select seq_region_id from seq_region where name = %s' %(chr))
        n = cursor.execute('select seq_region_id from %s.seq_region where name = %s' %(seq_region.db, chromosome))
        #print "Number of rows selected: ", n
        t = cursor.fetchall()[0]
        seq_region_ID = t[0]
        #print "Selected seq_region: ", seq_region_ID
        chr_seq = srdb[seq_region_ID]
        #print len(chr_seq)
        # convert an ensembl start coordinate to a Python zero-off coordinate
        python_start = start - 1
        # obtain the required interval
        ival = chr_seq[python_start:end]
        # If the required interval is on the reverse strand, then return the reverse and complimented sequence interval.
        if strand == -1:
            ival = -ival
        return ival

    def _fetch_units_by_seqregion(self, chromosome, start, end, strand, unit_name):
        '''Return a list of units (specified by unit_name, such as 'gene', 'transcript' or 'exon') found in a genomic region (defined by chromosome, start, end, strand).
        '''
        
        # create a SeqregionAdaptor object
        seq_region = self.getAdaptor('seq_region')
        # create a SeqregionAdaptor table object
        seq_regiontb = seq_region.tbobj
        
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        # create a SeqRegion object 
        srdb = SeqRegion(seq_regiontb, {17:hg18}, {17:'chr'})
        # create a [Exon/Gene/Transcript]Adaptor object depending on the specified unit_name
        unit_adaptor = self.getAdaptor(unit_name)
        # create a []Adaptor table object
        unit_TB = unit_adaptor.tbobj
        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(unit_TB, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))

        mapper = EnsemblMapper(annoDB, srdb) # find units for any sequence slice
        # retrieve the required seq_region_id from the seq_region table by chromosome name
        cursor = seq_region.cursor
        n = cursor.execute('select seq_region_id from %s.seq_region where name = %s' %(seq_region.db, chromosome))
        t = cursor.fetchall()[0]
        seq_region_ID = t[0]
        # obtain the chromosome sequence object based on the seq_region_ID
        chr_seq = srdb[seq_region_ID]
        # convert an ensembl start 1-offset coordinate to a Python zero-offset coordinate
        python_start = start - 1
        ival = chr_seq[python_start:end]    
        # find units in this sequence interval
        if strand == -1:
            unit_list = mapper[-ival]
        else:
            unit_list = mapper[ival]
        # test
        print '\n', unit_name, 'in interval', start, '-', end, 'on both strands of chromosome', chromosome,':'
        if len(unit_list) == 0:
            print 'None'
        else:
            print unit_list # find units in this interval on both strands
        
        return unit_list


class Adaptor(object):
    '''A base class that provides access to a generic table in a core ensembl   
    database'''
   
    def __init__(self, dbname, tbname, RowObj, cursor):
        self.db = dbname
        self.tb = tbname
        self.row = RowObj
        self.tbobj =  sqlgraph.SQLTable(self.db+'.'+self.tb, cursor, 
                                        itemClass=self.row)
        self.cursor = self.tbobj.cursor

    def __getitem__(self, i):
        return self.tbobj[i]


'''
class ExonTranscriptAdaptor(Adaptor):
    'Provides access to the exon_transcript table in an ensembl core database'

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'exon_transcript', sqlgraph.TupleO, cursor)
'''


class TranslationAdaptor(Adaptor):
    '''Provides access to the translation table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'translation', sqlgraph.TupleO, cursor)


class GeneStableIdAdaptor(Adaptor):
    '''Provides access to the gene_stable_id table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'gene_stable_id', sqlgraph.TupleO, cursor)


class TranscriptStableIdAdaptor(Adaptor):
    '''Provides access to the transcript_stable_id table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'transcript_stable_id', sqlgraph.TupleO, cursor)



class TranslationStableIdAdaptor(Adaptor):
    '''Provides access to the translation_stable_id table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'translation_stable_id', sqlgraph.TupleO, cursor)



class ExonStableIdAdaptor(Adaptor):
    '''Provides access to the exon_stable_id table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'exon_stable_id', sqlgraph.TupleO, cursor)



class SeqregionAdaptor(Adaptor):
    '''Provides access to the seq_region table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'seq_region', sqlgraph.TupleO, cursor)

    #def fetch_seqregion_by_something(self, something):

class ExonAdaptor(Adaptor):
    '''Provides access to the exon table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'exon', EnsemblRow, cursor)

    def fetch_exons_by_seqregion(self, chr, start, end, strand, driver):
        'find all the exons in a genomic region (defined by chromosoem, start, end, and strand)'
        
        # find all the exons on both strands of a genomic region
        unit_list =driver._fetch_units_by_seqregion(chr, start, end, strand, 'exon')
        # only return the list of exons found on the specified strand
        exons = []
        for annobj in unit_list:
            if annobj.orientation == 1:
                e = Exon(annobj.id)
                exons.append(e)
        if len(exons) == 0:
            print 'no exons found on this strand(', strand, ')'
        return exons

    def fetch_exons_by_transcriptID(self, transcript_id):
	cursor = self.cursor
        n = cursor.execute('select exon_id from %s.exon_transcript where transcript_id = %s' %(self.db, transcript_id))
        t = cursor.fetchall()
        exons = []
        if n == 0:
            return exons
        for row in t:
            e = Exon(row[0])
            exons.append(e)
        return exons

    def fetch_exons_by_translation(self, transcript_id, start_exon_id, end_exon_id):
        
        cursor = self.cursor
        n = cursor.execute('select rank from %s.exon_transcript where transcript_id = %%s and exon_id = %%s' %(self.db), (transcript_id, start_exon_id))
        t = cursor.fetchall()
        if n != 1:
            raise KeyError('Warning: duplicated!')
        start_rank = t[0][0]
        n = cursor.execute('select rank from %s.exon_transcript where transcript_id = %%s and exon_id = %%s' %(self.db), (transcript_id, end_exon_id))
        t = cursor.fetchall()
        if n != 1:
            raise KeyError('Warning: duplicated!')
        end_rank = t[0][0]
        
        n = cursor.execute('select exon_id from %s.exon_transcript where transcript_id = %%s and rank >= %%s and rank <= %%s' %(self.db), (transcript_id, start_rank, end_rank))
        t = cursor.fetchall()
        exons = []
        if n == 0:
            return exons
        for row in t:
            e = Exon(row[0])
            exons.append(e)
        return exons

class TranscriptAdaptor(Adaptor):
    'Provides access to the transcript table in an ensembl core database'

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'transcript', EnsemblRow, cursor)

    def fetch_transcripts_by_seqregion(self, chr, start, end, strand, driver):
        'find all the transcripts in a genomic region'

        # find all the transcripts on both strands of a genomic region
        unit_list = driver._fetch_units_by_seqregion(chr, start, end, strand, 'transcript')
        # return only the list of transcripts found on the specified strand
        transcripts = []
        for annobj in unit_list:
            if annobj.orientation == 1:
                t = Transcript(annobj.id)
                transcripts.append(t)
        if len(transcripts) == 0:
            print 'no transcript found on this strand(', strand, ')'
        return transcripts
    
    def fetch_transcripts_by_geneID(self, gene_id):
        'obtain all the transcripts that share the given gene_id'

        cursor = self.cursor
        #print 'select transcript_id from %s.transcript where gene_id = %s' %(self.db, gene_id)
        n = cursor.execute('select transcript_id from %s.transcript where gene_id = %%s' %(self.db), (gene_id))
        #n = cursor.execute('select transcript_id from %s.transcript where gene_id = %s' %(self.db, gene_id))

        t = cursor.fetchall()
        #print t
        transcripts = []
        if n == 0:
            return transcripts
        else:
            for row in t:
                t = Transcript(row[0])
                transcripts.append(t)
            return transcripts



class GeneAdaptor(Adaptor):
    '''Provides access to the gene table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'gene', EnsemblRow, cursor)

    def fetch_genes_by_seqregion(self, chr, start, end, strand, driver):
        'find all the genes in a genomic region'

        # find all the genes on both strands of a genomic region
        unit_list =driver._fetch_units_by_seqregion(chr, start, end, strand, 'gene')
        # return only the list of genes found on the specified strand
        genes = []
        for annobj in unit_list:
            if annobj.orientation == 1:
                g = Gene(annobj.id)
                genes.append(g)
        if len(genes) == 0:
            print 'no gene found on this strand(', strand, ')'
        return genes
    
    #def fetch_genes_by_something(self, something):
    
   
    
if __name__ == '__main__': # example code
    
    driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
    '''
    exon_adaptor = driver.getAdaptor('exon')
    #exons = exon_adaptor.fetch_exons_by_seqregion(1, 6023217, 6023986, 1, driver)
    #exons = exon_adaptor.fetch_exons_by_seqregion(10, 444866, 444957, -1, driver)
    exons = exon_adaptor.fetch_exons_by_seqregion(1, 1, 100000, 1, driver)
    for index, e in enumerate(exons):
       print '\nexon', index 
       e.getAttributes()
       
    gene_adaptor = driver.getAdaptor('gene')
    genes = gene_adaptor.fetch_genes_by_seqregion(1, 4274, 19669, -1, driver)
    for index, g in enumerate(genes):
        print '\ngene', index
        g.getAttributes()
        g.getSequence('gene')
    '''
    transcript_adaptor = driver.getAdaptor('transcript')
    '''
    transcripts = transcript_adaptor.fetch_transcripts_by_seqregion(1, 4274, 19669, -1, driver) 
    for index, t in enumerate(transcripts):
        print '\ntranscript', index
        t.getAttributes()
        #t.getSequence('transcript')
    '''
    print '\ntranscript_adaptor.fetch_transcripts_by_geneID(gene_id):'
    transcripts = transcript_adaptor.fetch_transcripts_by_geneID(34)
    if len(transcripts) == 0:
        print '\nNo transcript identified for this gene.'
    else:
        for index, t in enumerate(transcripts):
            print '\ntranscript ', index, ':'
            t.getAttributes()
            print 'length: ', len(t.getSequence('transcript'))
    
    '''
    #s = driver.fetch_sequence_by_region(1, 4274, 19669, 1)
    s = driver.fetch_sequence_by_region(10, 444866, 444957, 1)
    print "\nLength of the sequence: ", len(s)
    print "The sequence: ", str(s)
    '''


