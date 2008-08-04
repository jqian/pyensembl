from ensembl.seqregion import *
from ensembl.datamodel import *
import MySQLdb

class TranslationAdaptor(sqlgraph.SQLTable):
    '''Provides access to the translation table in an ensembl core database
    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> translationAdaptor = coreDBAdaptor.get_adaptor('translation')
    >>> translation = translationAdaptor[10]
    >>> translation.getAttributes()
    start_exon_id  =  61
    translation_id  =  10
    seq_start  =  98
    seq_end  =  767
    transcript_id  =  10
    end_exon_id  =  62
    id  =  10
    >>> translation.transcript_id
    10L
    >>> translation.start_exon_id
    61L
    >>> translation.end_exon_id
    62L
    '''
    def fetch_by_transcriptID(self, transcript_id):
        
        # retrieve row objects from the translation table, which statisfy the SQL select where clause
        t = self.tbobj.select('where transcript_id = %s', (transcript_id), None, 't1.translation_id') 
        translations = []
        for row in t:
            translation = Translation(row.translation_id)
            translations.append(translation)
        return translations



def get_registry(**kwargs):
    '''Provides a method to generate a connection to the ensembl server

    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> exonAdaptor = coreDBAdaptor.get_adaptor('exon')
    >>> exon = exonAdaptor[10]
    >>> exon.getAttributes()
    seq_region_start  =  5767
    seq_region_end  =  5810
    exon_id  =  10
    seq_region_id  =  226034
    is_current  =  1
    end_phase  =  0
    phase  =  1
    id  =  10
    seq_region_strand  =  -1
    '''

    if Registry._instance == None:
        Registry._instance = Registry(**kwargs)
    return Registry._instance

def _get_resource(resource_id):
    '''Obtain the required ensembl table from pygr.Data, return None if unable to find it in PYGRDATAPATH'''

    import pygr.Data
    try:
        return pygr.Data.getResource(resource_id)
    except pygr.Data.PygrDataNotFoundError:
        return None

def _save_resource(resource_id, tbobj):
    'Save the ensembl database table to pygr.Data'

    import pygr.Data
    tbobj.__doc__ = 'ensembl ' + resource_id.split('.')[3] + ' ' + resource_id.split('.')[4] + ' table'
    print tbobj.__doc__
    pygr.Data.addResource(resource_id, tbobj)
    # Save all pending data and schema to resource database
    pygr.Data.save()


class Registry(object):
    'Provide a connection to the ensembl server'

    _instance = None

    def __init__(self, **kwargs):
        self.kwargs = kwargs # connection arguments
        self.conn = sqlgraph.DBServerInfo(**self.kwargs)
        self.db_adaptor = {'core': CoreDBAdaptor}

    def get_DBAdaptor(self, db_species, db_type, db_version):
        'Obtain an adaptor to the given type of database'

        db_adaptor = self.db_adaptor[db_type]
        return db_adaptor(self, db_species, db_version)



"""
#class Adaptor(object):
class Adaptor(sqlgraph.SQLTable):
    '''A base class that provides access to a generic table in a core ensembl   
    database'''
   
    def __init__(self, dbname, tbname, RowObj, conn):
        
        sqlgraph.SQLTable.__init__(self, dbname+'.'+tbname, serverInfo=conn, itemClass=RowObj)

        # Get the tbobj from pygr.Data
        self.resource_id = 'Bio.MySQL.EnsemblTB.' + self.db + '.' + self.tb
        #print self.resource_id
        self.tbobj = _get_Resource(self.resource_id)
        if self.tbobj == None:
            # Create a pickleable SQLTable object
            self.tbobj = sqlgraph.SQLTable(self.db+'.'+self.tb, serverInfo=conn,itemClass=self.row)
            #self.tbobj = sqlgraph.SQLTable(self.db+'.'+self.tb, serverInfo=conn)
            # Save self.tbobj to pygr.Data
            #_save_resource(self.resource_id, self.tbobj)
        self.cursor = self.tbobj.cursor

    def __getitem__(self, i):
        'Create an object based on an item from an ensembl database table'

        item = self.tbobj[i]
        itemClass = self.row_object[self.tb]
        itemObj = itemClass(item)
        return itemObj
"""   
        
'''
class ExonTranscriptAdaptor(Adaptor):
    'Provides access to the exon_transcript table in an ensembl core database'

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'exon_transcript', sqlgraph.TupleO, cursor)
'''


class MetaCoordAdaptor(sqlgraph.SQLTable):
    '''Provides access to the meta_coord table in an ensembl core database'''
    

class XrefAdaptor(sqlgraph.SQLTable):
    '''Provides access to the xref table in an ensembl core database'''




class GeneStableIdAdaptor(sqlgraph.SQLTable):
    '''Provides access to the gene_stable_id table in an ensembl core database'''

   

    def fetch_by_stable_id(self, stable_id):
        t = self.tbobj.select('where stable_id = %s', (stable_id), None, 't1.gene_id')
        genes = []
        for row in t:
            gene_id = row.gene_id
            gene = Gene(gene_id)
            genes.append(gene)
        return genes


class TranscriptStableIdAdaptor(sqlgraph.SQLTable):
    '''Provides access to the transcript_stable_id table in an ensembl core database'''


    def fetch_by_stable_id(self, stable_id):
        t = self.tbobj.select('where stable_id = %s', (stable_id), None, 't1.transcript_id')
        transcripts = []
        for row in t:
            transcript_id = row.transcript_id
            transcript = Transcript(transcript_id)
            transcripts.append(transcript)
        return transcripts



class TranslationStableIdAdaptor(sqlgraph.SQLTable):
    '''Provides access to the translation_stable_id table in an ensembl core database'''

    
    def fetch_by_stable_id(self, stable_id):
        t = self.tbobj.select('where stable_id = %s', (stable_id), None, 't1.translation_id')
        translations = []
        for row in t:
            translation_id = row.translation_id
            translation = Translation(translation_id)
            translations.append(translation)
        return translations




class ExonStableIdAdaptor(sqlgraph.SQLTable):
    '''Provides access to the exon_stable_id table in an ensembl core database'''


    def fetch_by_stable_id(self, stable_id):
        t = self.tbobj.select('where stable_id = %s', (stable_id), None, 't1.exon_id')
        exons = []
        for row in t:
            exon_id = row.exon_id
            exon = Exon(exon_id)
            exons.append(exon)
        return exons



class SeqregionAdaptor(sqlgraph.SQLTable):
    '''Provides access to the seq_region table in an ensembl core database'''

    

class PeptideArchiveAdaptor(sqlgraph.SQLTable):
    '''Provides access to the peptide_archive table in an ensembl core database'''

   

"""
class GeneArchiveAdaptor(Adaptor):
    '''Provides access to the gene_archive table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'gene_archive', sqlgraph.TupleO, cursor)
"""


class PredictionExonAdaptor(sqlgraph.SQLTable):
    '''Provides access to the prediction_exon table in an ensembl core database'''

   
    def fetch_all_by_ptranscriptID(self, predic_transcriptID):
        t = self.tbobj.select('where prediction_transcript_id = %s', (predic_transcriptID)) 
        predic_exons = []
        for row in t:
            predic_exon = PredictionExon(row)
            predic_exons.append(predic_exon)
        return predic_exons


class PredictionTranscriptAdaptor(sqlgraph.SQLTable):
    '''Provides access to the prediction_transcript table in an ensembl core database'''


    def fetch_by_display_label(self, display_label):
        '''Obtain prediction_transcript objects by a given display_label
 
        >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> ptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
        >>> pts = ptAdaptor.fetch_by_display_label('GENSCAN00000036948')
        >>> for index, pt in enumerate(pts):
        ...     print 'prediction_transcript', index, ':'
        ...     pt.getAttributes()
        ...     print 'start:', pt.start
        ...     print 'seq_region_start:', pt.seq_region_start
        ...
        prediction_transcript 0 :
        analysis_id  =  5
        seq_region_start  =  48849
        seq_region_end  =  66632
        seq_region_id  =  149878
        prediction_transcript_id  =  36948
        display_label  =  GENSCAN00000036948
        id  =  36948
        seq_region_strand  =  -1
        start: 48848
        seq_region_start: 48849

        '''
    
        #t = self.tbobj.select('where display_label = %s', (display_label), EnsemblRow)
        t = self.select('where display_label = %s', (display_label))
        prediction_transcripts = []
        for row in t:
            prediction_transcript = self[row.id]
            prediction_transcripts.append(prediction_transcript)
        return prediction_transcripts
            
    def fetch_all_by_seqregion(self, chr, start, end, strand, driver):
        'find all the prediction_transcripts in a genomic region'

        # find all the prediction_transcripts on both strands of a genomic region
        unit_list = driver._fetch_units_by_seqregion(chr, start, end, strand, 'transcript')
        # return only the list of transcripts found on the specified strand
        prediction_transcripts = []
        for annobj in unit_list:
            if annobj.orientation == 1:
                t = PredictionTranscript(annobj.id)
                prediction_transcripts.append(t)
        if len(prediction_transcripts) == 0:
            print 'no prediction transcript found on this strand(', strand, ')'
        return prediction_transcripts

    def fetch_all_by_seqregion(self, chr, start, end, strand, driver):
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
    

class ExonAdaptor(sqlgraph.SQLTable):
    '''Provides access to the exon table in an ensembl core database'''

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


class TranscriptAdaptor(sqlgraph.SQLTable):
    'Provides access to the transcript table in an ensembl core database'

   
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



class GeneAdaptor(sqlgraph.SQLTable):
    '''Provides access to the gene table in an ensembl core database'''


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

class CoreDBAdaptor(object):
    'Provide access to a core database in the ensembl database system'

    # static attributes
    dbType = 'core'
    TBAdaptorClass = {'exon': ExonAdaptor, 'gene': GeneAdaptor, 'transcript': TranscriptAdaptor, 'seq_region': SeqregionAdaptor, 'translation': TranslationAdaptor, 'gene_stable_id': GeneStableIdAdaptor, 'transcript_stable_id': TranscriptStableIdAdaptor, 'translation_stable_id': TranslationStableIdAdaptor, 'exon_stable_id': ExonStableIdAdaptor, 'xref': XrefAdaptor, 'prediction_transcript': PredictionTranscriptAdaptor, 'meta_coord': MetaCoordAdaptor, 'prediction_exon': PredictionExonAdaptor, 'peptide_archive': PeptideArchiveAdaptor}
    #TBAdaptorClass = {'translation': TranslationAdaptor}
    RowClass = {'exon': Exon, 'gene': Gene, 'transcript': Transcript, 'seq_region': Seqregion, 'translation': Translation, 'gene_stable_id': GeneStableID, 'transcript_stable_id': TranscriptStableID, 'translation_stable_id': TranslationStableID, 'exon_stable_id': ExonStableID, 'xref': Xref, 'prediction_transcript': PredictionTranscript, 'prediction_exon': PredictionExon, 'peptide_archive': PeptideArchive}
    def __init__(self, registry, dbSpecies, dbVersion):
        
        # instance attributes
        self.conn = registry.conn
        self.dbSpecies = dbSpecies
        self.dbVersion = dbVersion
    
    def get_adaptor(self, tbname): 
        'Obtain a particular table adaptor of a core database table'

        # Get the tbobj from pygr.Data
        dbname = self.dbSpecies+'_'+CoreDBAdaptor.dbType+'_'+self.dbVersion+'.'+tbname
        resource_id = 'Bio.MySQL.EnsemblTB.' + dbname
        #print self.resource_id
        tbAdaptor = _get_resource(resource_id)
        if tbAdaptor == None:
            # Create a pickleable SQLTable adaptor
            adaptorClass = CoreDBAdaptor.TBAdaptorClass[tbname]
            rowClass = CoreDBAdaptor.RowClass[tbname]
            #self.tbobj = sqlgraph.SQLTable(self.db+'.'+self.tb, serverInfo=conn)
            #self.tbobj = sqlgraph.SQLTable(self.db+'.'+self.tb, serverInfo=conn,itemClass=self.row)
            tbAdaptor = adaptorClass(dbname, itemClass=rowClass, serverInfo=self.conn)
           
            # Save self.tbobj to pygr.Data
            #_save_resource(resource_id, tbAdaptor)
        #self.cursor = self.tbobj.cursor
        return tbAdaptor
        #return adaptor_name(self.db_species + '_' + self.db_type + '_' + self.db_version, self.conn)

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
        # t = seq_regiontb.select('where name = %s', (chromosome), None, 't1.seq_region)
        t = seq_regiontb.select('where name = %s', (chromosome))
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
        
        # human genome
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() 

        # ensembl sequences
        tbname = self.dbname + '.dna'
        dna = sqlgraph.SQLTable(tbname, cursor, 
                            itemClass=EnsemblDNA,
                            itemSliceClass=seqdb.SeqDBSlice,
                            attrAlias=dict(seq='sequence'))


        # create a SeqRegion object 
        srdb = SeqRegion(seq_regiontb, {17:hg18, 4:dna}, {17:'chr', 4:None})
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

        # retrieve the slice according to the specified strand
        if strand == -1:
            slice = -ival
        else:
            slice = ival

        # check which coord_system the required feature belongs to
        meta_coord_adaptor = self.getAdaptor('meta_coord')
        meta_coord_tb = meta_coord_adaptor.tbobj
        t = meta_coord_tb.select('where table_name = %s', (unit_name))
        row = t.next()
        unit_cs_id = row.coord_system_id
        if unit_cs_id == 4:
            # transform the sequence interval from the genomic coordinates to the contig coordinates
        
            # create an assembly mapper
            amap = AssemblyMapper(srdb, 4, 17)         
            # map the genomic slice to a contig slice
            contig_slice = (~amap)[slice]         
            #contig = srdb[149878][pt.sequence.start:pt.sequence.stop]
            slice = contig_slice
        
        # find units in this genomic or contig sequence interval
        unit_list = mapper[slice]

        # test
        print '\n', unit_name, 'in interval', start, '-', end, 'on both strands of chromosome', chromosome,':'
        if len(unit_list) == 0:
            print 'None'
        else:
            print unit_list # find units in this interval on both strands

       
        return unit_list


            
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
    
def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__': # example code
    
    _test()
    '''
    serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    '''
    '''
    translationAdaptor = coreDBAdaptor.get_adaptor('translation')
    translation = translationAdaptor[10]
    translation.getAttributes()
    #exons = translation.getExons()
    #print 'id', 'start'
    #for e in exons:
    #    print e.id, e.start
    '''
    '''
    ptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
    pts = ptAdaptor.fetch_by_display_label('GENSCAN00000036948')
    
    for index, pt in enumerate(pts):
        print 'prediction_transcript', index, ':'
        pt.getAttributes()
        print 'start:', pt.start
        print 'seq_region_id:', pt.seq_region_id
        print 'seq_region_start:', pt.seq_region_start
        print 'seq_region_end:', pt.seq_region_end
        print 'seq_region_strand:', pt.seq_region_strand
        print 'analysis_id:', pt.analysis_id
        print 'display_label:', pt.display_label
    '''
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
    
    
    prediction_transcript_adaptor = driver.getAdaptor('prediction_transcript')
    prediction_transcripts = prediction_transcript_adaptor.fetch_by_display_label('GENSCAN00000036948')
    for index, pt in enumerate(prediction_transcripts):
        print 'prediction_transcript', index, ':'
        pt.getAttributes()
        print 'start:', pt.rowobj.start
        print 'seq_region_id:', pt.getSeqregionID()
        print 'seq_region_start:', pt.getSeqregionStart()
        print 'seq_region_end:', pt.getSeqregionEnd()
        print 'seq_region_strand:', pt.getOrientation()
        print 'analysis_id:', pt.getAnalysisID()
        print 'display_label:', pt.getDisplayLabel()
    '''   

