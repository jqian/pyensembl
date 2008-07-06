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
        self.tb_adaptor = {'exon': ExonAdaptor, 'gene': GeneAdaptor, 'transcript': TranscriptAdaptor, 'seq_region': SeqregionAdaptor, 'translation': TranslationAdaptor, 'gene_stable_id': GeneStableIdAdaptor, 'transcript_stable_id': TranscriptStableIdAdaptor, 'translation_stable_id': TranslationStableIdAdaptor, 'exon_stable_id': ExonStableIdAdaptor, 'xref': XrefAdaptor}

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
        #print srdb
        # retrieve the required seq_region_id from the seq_region table by chromosome name
        '''
        cursor = seq_region.cursor
        #n = srdb.seqRegionDB.cursor.execute('select seq_region_id from seq_region where name = %s' %(chr))
        
        n = cursor.execute('select seq_region_id from %s.seq_region where name = %s' %(seq_region.db, chromosome))
        #print "Number of rows selected: ", n
        t = cursor.fetchall()[0]
        seq_region_ID = t[0]
        '''
        t = seq_regiontb.select('where name = %s', (chromosome), None, 't1.seq_region_id')
        # t = seq_regiontb.select('where name = %s', (chromosome))
        #seq_region_ID = t.next().seq_region_id
        row = t.next()
        seq_region_ID = row.seq_region_id
        #print srdb
        
        #print "Selected seq_region: ", seq_region_ID
        #print seq_regiontb[seq_region_ID].coord_system_id
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
        '''
        cursor = seq_region.cursor
        n = cursor.execute('select seq_region_id from %s.seq_region where name = %s' %(seq_region.db, chromosome))
        t = cursor.fetchall()[0]
        seq_region_ID = t[0]
        '''
        t = seq_regiontb.select('where name = %s', (chromosome), None, 't1.seq_region_id')
        row = t.next()
        seq_region_ID = row.seq_region_id
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


class XrefAdaptor(Adaptor):
    '''Provides access to the xref table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'xref', sqlgraph.TupleO, cursor)


class TranslationAdaptor(Adaptor):
    '''Provides access to the translation table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'translation', sqlgraph.TupleO, cursor)

    def fetch_translations_by_transcriptID(self, transcript_id):
        
        # retrieve row objects from the translation table, which statisfy the SQL select where clause
        t = self.tbobj.select('where transcript_id = %s', (transcript_id), None, 't1.translation_id') 
        translations = []
        for row in t:
            translation = Translation(row.translation_id)
            translations.append(translation)
        return translations



class GeneStableIdAdaptor(Adaptor):
    '''Provides access to the gene_stable_id table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'gene_stable_id', sqlgraph.TupleO, cursor)

    def fetch_by_stable_id(self, stable_id):
        t = self.tbobj.select('where stable_id = %s', (stable_id), None, 't1.gene_id')
        genes = []
        for row in t:
            gene_id = row.gene_id
            gene = Gene(gene_id)
            genes.append(gene)
        return genes


class TranscriptStableIdAdaptor(Adaptor):
    '''Provides access to the transcript_stable_id table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'transcript_stable_id', sqlgraph.TupleO, cursor)

    def fetch_by_stable_id(self, stable_id):
        t = self.tbobj.select('where stable_id = %s', (stable_id), None, 't1.transcript_id')
        transcripts = []
        for row in t:
            transcript_id = row.transcript_id
            transcript = Transcript(transcript_id)
            transcripts.append(transcript)
        return transcripts



class TranslationStableIdAdaptor(Adaptor):
    '''Provides access to the translation_stable_id table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'translation_stable_id', sqlgraph.TupleO, cursor)
    
    def fetch_by_stable_id(self, stable_id):
        t = self.tbobj.select('where stable_id = %s', (stable_id), None, 't1.translation_id')
        translations = []
        for row in t:
            translation_id = row.translation_id
            translation = Translation(translation_id)
            translations.append(translation)
        return translations




class ExonStableIdAdaptor(Adaptor):
    '''Provides access to the exon_stable_id table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'exon_stable_id', sqlgraph.TupleO, cursor)

    def fetch_by_stable_id(self, stable_id):
        t = self.tbobj.select('where stable_id = %s', (stable_id), None, 't1.exon_id')
        exons = []
        for row in t:
            exon_id = row.exon_id
            exon = Exon(exon_id)
            exons.append(exon)
        return exons



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
        
        # Obtain row objects that satisfy the select where clause
        t = self.tbobj.select('where gene_id = %s', (gene_id), None, 't1.transcript_id')
        transcripts = []
        for row in t:
            t = Transcript(row.transcript_id)
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
    
    def fetch_genes_by_externalRef(self, external_ref_label):
        '''Return all the genes that are associated with the given external reference or an empty set.'''
        
        '''
        cursor = self.cursor
        #print('select xref_id from %s.xref where display_label = %s' %(self.db, external_ref_label)) 
        
        n = cursor.execute('select xref_id from %s.xref where display_label = %%s' %(self.db), (external_ref_label))
        t = cursor.fetchall()
        '''
        # Retrieve xref_id(s) that are associated with the given external reference
        
        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        xref_adaptor = driver.getAdaptor('xref')
        t = xref_adaptor.tbobj.select('where display_label = %s', (external_ref_label), None, 't1.xref_id')
        
        # Retrieve gene_id(s) that are associated with the returned xref_id(s)
        geneIDs = []
        for row in t:
            xref_id = row.xref_id
            #print xref_id
            gene_rows = self.tbobj.select('where display_xref_id = %s', (xref_id), None, 't1.gene_id')
            for gid_row in gene_rows:
                geneIDs.append(gid_row.gene_id)
        # Create genes based on the retrieved gene_ids
        genes = []
        for gid in geneIDs:
            g = Gene(gid)
            genes.append(g)

        return genes
            
def _retrieve_units_tester(units, unit_name):
    for index, u in enumerate(units):
        print unit_name, index, ':'
        u.getAttributes()
        s = u.getSequence(unit_name)
        print '\nLength:', len(s)
        print '\nSequence:', str(s)

def _fetch_by_stableID_tester(myAdaptor, stable_id):
    units = myAdaptor.fetch_by_stable_id(stable_id)
    for index, u in enumerate(units):
        print index, ':'
        u.getAttributes()
    
   
    
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
    
    print '\ngene_adaptor.fetch_genes_by_externalRef():'
    gene_adaptor = driver.getAdaptor('gene')
    genes = gene_adaptor.fetch_genes_by_externalRef('IQSEC3')
    #genes = gene_adaptor.fetch_genes_by_externalRef('GO:0000228')
    if len(genes) == 0:
        print '\nNo gene is associated with the given external reference label.'
    else:
        for index, g in enumerate(genes):
            print '\ngene', index, ':'
            g.getAttributes()
       
    genes = gene_adaptor.fetch_genes_by_seqregion(1, 4274, 19669, -1, driver)
    for index, g in enumerate(genes):
        print '\ngene', index
        g.getAttributes()
        g.getSequence('gene')
    
    transcript_adaptor = driver.getAdaptor('transcript')
    
    transcripts = transcript_adaptor.fetch_transcripts_by_seqregion(1, 4274, 19669, -1, driver) 
    for index, t in enumerate(transcripts):
        print '\ntranscript', index
        t.getAttributes()
        #t.getSequence('transcript')
    
    print '\ntranscript_adaptor.fetch_transcripts_by_geneID(gene_id):'
    transcripts = transcript_adaptor.fetch_transcripts_by_geneID(8946)
    if len(transcripts) == 0:
        print '\nNo transcript identified for this gene.'
    else:
        for index, t in enumerate(transcripts):
            print '\ntranscript ', index, ':'
            t.getAttributes()
            print 'length: ', len(t.getSequence('transcript'))
    
    
    #s = driver.fetch_sequence_by_region(1, 4274, 19669, 1)
    s = driver.fetch_sequence_by_region(10, 444866, 444957, 1)
    print "\nLength of the sequence: ", len(s)
    print "The sequence: ", str(s)
   
    '''
    gene_stableID_adaptor = driver.getAdaptor('gene_stable_id')
    print "\ntest gene_stableID_adaptor.fetch_by_stable_id('ENSG00000215911')"
    #genes = gene_stableID_adaptor.fetch_by_stable_id('ENSG00000215911')
    #_retrieve_units_tester(genes, 'gene')
    _fetch_by_stableID_tester(gene_stableID_adaptor, 'ENSG00000215911')

    transcript_stableID_adaptor = driver.getAdaptor('transcript_stable_id')
    print "\ntest transcript_stableID_adaptor.fetch_by_stable_id('ENST00000382841')"
    _fetch_by_stableID_tester(transcript_stableID_adaptor, 'ENST00000382841')

    translation_stableID_adaptor = driver.getAdaptor('translation_stable_id')
    print "\ntest translation_stableID_adaptor.fetch_by_stable_id('ENSP00000317958')"
    _fetch_by_stableID_tester(translation_stableID_adaptor, 'ENSP00000317958')

    exon_stableID_adaptor = driver.getAdaptor('exon_stable_id')
    print "\ntest exon_stableID_adaptor.fetch_by_stable_id('ENSE00001493538')"
    _fetch_by_stableID_tester(exon_stableID_adaptor, 'ENSE00001493538')


