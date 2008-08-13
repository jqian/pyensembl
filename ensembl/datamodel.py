from ensembl.adaptor import *
from ensembl.seqregion import *
from seqregion import EnsemblRow


class BaseModel(sqlgraph.TupleO):
    '''A generic interface to an item object in a table in the ensembl database
    '''
    
    '''    
    def getAttributes(self):
        'print out this row record'
	#for k, v in self.rowobj._attrcol.iteritems():
        #    print k, ' = ', self.rowobj.data[v]
        for k, v in self._attrcol.iteritems():
	    print k, ' = ', self.data[v]
    '''

    # test
    def getDBName(self):
        return self.db.name

class Dna(EnsemblDNA):
    '''An interface to a row in the dna table'''


class CoordSystem(BaseModel):
    '''An interface to a row in the coord_system table'''


class Xref(BaseModel):
    '''An interface to a record in the xref table in any ensembl core database'''

    
    



    


class GeneStableID(BaseModel):
    '''An interface to a record in the gene_stable_id table in any ensembl core database'''


class TranscriptStableID(BaseModel):
    '''An interface to a record in the transcript_stable_id table in any ensembl core database'''

class ExonStableID(BaseModel):
    '''An interface to a record in the exon_stable_id table in any ensembl core database'''

class TranslationStableID(BaseModel):
    '''An interface to a record in the translation_stable_id table in any ensembl core database'''

class PeptideArchive(BaseModel):
    '''
    An interface to a row record in the peptide_archive table
    
    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> pepAdaptor = coreDBAdaptor.get_adaptor('peptide_archive')
    >>> peptide = pepAdaptor[100]
    >>> print peptide.peptide_seq
    MKKARNDEYENLFNMIVEIPRWTNAKMEIATKEPMNPIKQYVKDGKLRYVANIFPYKGYIWNYGTLPQTWEDPHEKDKSTNCFGDNDPIDVCEIGSKILSCGEVIHVKILGILALIDEGETDWKLIAINANDPEASKFHDIDDVKKFKPGYLEATLNWFRLYKVPDGKPENQFAFNGEFKNKAFALEVIKSTHQCWKALLMKKCNGGAINCTNVQISDSPFRCTQEEARSLVESVSSSPNKESNEEEQVWHFLGK

    '''




class PredictionExon(BaseModel, EnsemblRow):
    '''
    An interface to a row record in the prediction_exon table

    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> pexonAdaptor = coreDBAdaptor.get_adaptor('prediction_exon')
    >>> pExon = pexonAdaptor[10]
    >>> print pExon.start
    31473
    >>> print pExon.seq_region_start
    31474
    >>> print pExon.prediction_transcript_id
    3
    >>> ptranscript = pExon.get_prediction_transcript()
    >>> print ptranscript.id
    3

    '''

    def get_prediction_transcript(self):
        'Obtain the prediction transcript the given prediction exon belongs to'
        
        ptranscriptPExons = _get_featureMapper('prediction_transcript', 'prediction_exon', self.db.name)
        ptranscript = (~ptranscriptPExons)[self]
        return ptranscript

    
   
class PredictionTranscript(BaseModel, EnsemblRow):
    '''
    An interface to a prediction_transcript record in any ensembl core database
    
    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> ptranscriptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
    >>> ptranscript = ptranscriptAdaptor[150]
    >>> print ptranscript.start
    18007
    >>> print ptranscript.seq_region_start
    18008
    >>> pexons = ptranscript.get_prediction_exons()
    >>> for pe in pexons:
    ...     print pe.id
    ...
    822
    823
    824
    '''

    def get_prediction_exons(self):
        'Obtain all the prediction exons for this prediction transcript'

        ptranscriptPExons = _get_featureMapper('prediction_transcript', 'prediction_exon', self.db.name)
        pexons = ptranscriptPExons[self]
        return pexons


class Seqregion(BaseModel):
    '''An interface to a seq_region record in the seq_region table in any 
    ensembl core database'''



class StableObj(BaseModel):
    '''An interface to a generic rowObj in a table that has an ensembl stable_idassigned to it.  Subclasses of this class are Gene, Transcript, Translation and Exon'''
    """
    def _getSeqregionDB(self):      
        'a private helper method to create a seq_region database object'

        seq_regiontb = self.driver.getAdaptor('seq_region').tbobj
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        srdb = SeqRegion(seq_regiontb, {17:hg18}, {17:'chr'})
        return srdb 


    def _getAnnotationDB(self, unit_tbname):
        'a private helper method to create an annotation db object'

        unit_TB = self.driver.getAdaptor(unit_tbname).tbobj
        # obtain a seq_region db object
        srdb = self._getSeqregionDB()
        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(unit_TB, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))
        return annoDB

    def _getSeqregionSeq(self):
        '''obtain the sequence of a seq_region for a sliceable object.
        This seq_region is defined by seq_region_id, seq_region_start, seq_region_end and seq_region_strand.
        '''
        
        srdb = self._getSeqregionDB()
        seq_region_ID = self.getSeqregionID()
        aSeqregion = srdb[seq_region_ID]
        aStart = self.rowobj.start 
        aEnd = self.getSeqregionEnd()
        aStrand = self.getOrientation()
        if aStrand == -1:
            aSequence = -aSeqregion[aStart:aEnd]
        else:
            aSequence = aSeqregion[aStart:aEnd]
        return aSequence
     
    def getSequence(self, unit_tbname):
        '''Get a sequence object of a sliceable object.
        Note:  the sequence object is retrieved from the strand it actually resides on
        '''
        
        # obtain an annotation database of a sliceable table object
        annoDB = self._getAnnotationDB(unit_tbname)
        
        unitID = self.rowobj.id
        unitobj = annoDB[unitID]
        s = unitobj.sequence
        return s
    """
    def _get_stable_id_obj(self):
        'Get a rowobj from a particular stable_id table'

        # get the database adaptor 
        name = self.db.name
        dbName = name.split('.')[0]
        dbSpecies = dbName.split('_')[0] + '_' + dbName.split('_')[1]
        dbType = dbName.split('_')[2]
        dbVersion = dbName.split('_')[3] + '_' + dbName.split('_')[4]
        from ensembl.adaptor import _get_DB_adaptor
        dbAdaptor = _get_DB_adaptor(dbSpecies, dbType, dbVersion)
        
        # get the particular stable_id table adaptor
        tbName = name.split('.')[1]
        stableIDtbName = tbName + '_' + 'stable_id'
        stableIDAdaptor = dbAdaptor.get_adaptor(stableIDtbName)
        sliceableID = self.id
        stableIDObj = stableIDAdaptor[sliceableID]
        return stableIDObj 
        
    def get_stable_id(self):
        '''Return the stable_id of this Sliceable, none if it doesn't have one'''
        stableIDObj = self._get_stable_id_obj()
        if stableIDObj == None:
            return None
        else:
            stableID = stableIDObj.stable_id
            return stableID

    def get_created_date(self):
        'Obtain the ensembl created date for the Sliceable if it has a stable id'
        stableIDObj = self._get_stable_id_obj()
        if stableIDObj == None:
            return None
        else:
            createdDate = stableIDObj.created_date
            return createdDate
    
    def get_modified_date(self):
        'Obtain the ensembl modified date for the Sliceable if it has a stable id'
        
        stableIDObj = self._get_stable_id_obj()
        if stableIDObj == None:
            return None
        else:
            modifiedDate = stableIDObj.modified_date
            return modifiedDate

    def get_version(self):
        'Obtain the ensembl version for the Sliceable if it has a stable id'
        
        stableIDObj = self._get_stable_id_obj()
        if stableIDObj == None:
            return None
        else:
            version = stableIDObj.version
            return version

        

    

class Exon(StableObj, EnsemblRow):
    '''An interface to an exon record in the exon table in an ensembl core 
    database


    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> exonAdaptor = coreDBAdaptor.get_adaptor('exon')
    >>> exon = exonAdaptor[95160]
    >>> print exon.get_stable_id()
    ENSE00001493538
    >>> print exon.get_created_date()
    2006-03-10 00:00:00
    >>> print exon.get_modified_date()
    2006-03-10 00:00:00
    >>> print exon.get_version()
    1
    >>> transcripts = exon.get_all_transcripts()
    >>> for t in transcripts:
    ...     print t.id, '       ', t.get_stable_id()
    ...
    15960         ENST00000382841
    '''
    #def getGene(self):


    def get_all_transcripts(self):
        'Obtain all the transcripts this exon belongs to'

        
        transcriptExons = _get_featureMapper('transcript', 'exon', self.db.name)
        transcripts = (~transcriptExons)[self]

        return transcripts

    
    # test
    def get_start(self):
        return self.start

def _get_featureMapper(sourceDBName, targetDBName, name):
    
    dbName = name.split('.')[0]
    dbSpecies = dbName.split('_')[0] + '_' + dbName.split('_')[1]
    dbType = dbName.split('_')[2]
    dbVersion = dbName.split('_')[3] + '_' + dbName.split('_')[4]
    from ensembl.adaptor import _get_DB_adaptor
    dbAdaptor = _get_DB_adaptor(dbSpecies, dbType, dbVersion)
    mapper = dbAdaptor._fetch_featureMapper(sourceDBName, targetDBName)
    return mapper
    

class Transcript(StableObj, EnsemblRow):
    '''An interface to a transcript record in the transcript table in an ensembl core database
    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
    >>> transcript = transcriptAdaptor[15960]
    >>> print transcript.get_stable_id()
    ENST00000382841
    >>> print transcript.get_created_date()
    2006-03-10 00:00:00
    >>> print transcript.get_modified_date()
    2006-03-10 00:00:00
    >>> print transcript.get_version()
    1
    >>> exons = transcript.get_all_exons()
    >>> for e in exons:
    ...     print e.id, ' ', e.get_stable_id()
    ...
    95144   ENSE00001493540
    95152   ENSE00001493539
    95160   ENSE00001493538
    95020   ENSE00000893353
    95035   ENSE00000893354
    95050   ENSE00000893355
    95059   ENSE00000893356
    95069   ENSE00000893357
    95081   ENSE00000893358
    95088   ENSE00000893359
    95101   ENSE00000893360
    95110   ENSE00000893361
    95172   ENSE00001493527
    >>> gene = transcript.get_gene()
    >>> print gene.id, gene.get_stable_id()
    8946 ENSG00000120645
    >>> translation = transcript.get_translation()
    >>> print translation.id, translation.get_stable_id()
    15121 ENSP00000372292

    '''

    def get_all_exons(self):
        'obtain all the exons of this transcript'

        
        transcriptExons = _get_featureMapper('transcript', 'exon', self.db.name)
        exons = transcriptExons[self]

        return exons

    #def getExternalRefs(self):
    #    return getExternalRefs(true);

    #def getExternalRefs(self, includeTranslation):

    def get_gene(self):
        'obtain its gene'
        
        geneTranscripts = _get_featureMapper('gene', 'transcript', self.db.name)
        gene = (~geneTranscripts)[self]
        return gene

    def get_translation(self):
        'obtain its translation'
        
        transcriptTranslation = _get_featureMapper('transcript', 'translation', self.db.name)
        translation = transcriptTranslation[self]
        return translation

    #def getSupportingFeatures(self):

    #def getAnalysis(self):


class Gene(StableObj, EnsemblRow):
    '''An interface to a gene record in the gene table in an ensembl core database
    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> geneAdaptor = coreDBAdaptor.get_adaptor('gene')
    >>> gene = geneAdaptor[8946]
    >>> transcripts = gene.get_transcripts()
    >>> for t in transcripts:
    ...     print t.id, t.get_stable_id()
    ...
    15950 ENST00000326261
    15960 ENST00000382841

    '''

       
    def get_transcripts(self):
        'return transcripts if available, otherwise empty list'
        
        geneTranscripts = _get_featureMapper('gene', 'transcript', self.db.name)
        transcripts = geneTranscripts[self]
        return transcripts
        
    def getExternalRefs(self):
        '''References to external databases.  They are a collection of all the 
        transcript.externalRefs.
        Returns a list of ExternalRef objects, list is empty if no external 
        references available.
        '''
        return getExternalRefs(self, true)

    def getExternalRefs(self, includeTranscriptsAndTranslations):
        'includeTranscriptsAndTranslations: a boolean flag'


class Translation(StableObj):
    '''An interface to an item in the translation table in any ensembl core database
    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> translationAdaptor = coreDBAdaptor.get_adaptor('translation')
    >>> translation = translationAdaptor[15121]
    >>> transcript = translation.get_transcript()
    >>> print transcript.id, transcript.get_stable_id()
    15960 ENST00000382841

    '''

    def get_transcript(self):
        'obtain its transcript'
        
        transcriptTranslation = _get_featureMapper('transcript', 'translation', self.db.name)
        transcript = (~transcriptTranslation)[self]
        return transcript

    def getExons(self):
        transcript_id = self.rowobj.transcript_id
        start_exon_id = self.get_start_exon_id()
        end_exon_id = self.get_end_exon_id()
        driver = self.driver
        exon_adaptor = driver.getAdaptor('exon')
        exons = exon_adaptor.fetch_exons_by_translation(transcript_id, start_exon_id, end_exon_id)
        return exons


def _getExons_tester(exons):
    'retrieve and print out all the exons returned by the *.getExons()'  
  
    if len(exons) == 0:
        print '\nno exon returned'
    else:
        for index, e in enumerate(exons):
            print '\nexon ', index, ':'
            e.getAttributes()
            #print 'Sequence:', str(e.getSequence())
            print 'Length of the sequence:', len(e.getSequence('exon'))

def _test():
    import doctest
    doctest.testmod()

if __name__ == '__main__': # example code
    
    _test()
    '''
    serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    
    pexonAdaptor = coreDBAdaptor.get_adaptor('prediction_exon')
    pExon = pexonAdaptor[10]
    print pExon.start
    print pExon.seq_region_start
    print pExon.prediction_transcript_id
    ptranscript = pExon.get_prediction_transcript()
    print ptranscript.id

    ptranscriptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
    ptranscript = ptranscriptAdaptor[3]
    print ptranscript.start
    print ptranscript.seq_region_start
    pexons = ptranscript.get_prediction_exons()
    for pe in pexons:
        print pe.id
    
    exonAdaptor = coreDBAdaptor.get_adaptor('exon')
    exon = exonAdaptor[95160]
    exon.get_stable_id()
    exon.get_created_date()
    exon.get_modified_date()
    exon.get_version()
    transcripts = exon.get_all_transcripts()
    print 'transcript_id', 'transcript_stable_id'
    for t in transcripts:
        print t.id, '       ', t.get_stable_id()
    
    transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
    transcript = transcriptAdaptor[15960]
    #transcript.get_stable_id()
    #transcript.get_created_date()
    #transcript.get_modified_date()
    #transcript.get_version()
    #exons = transcript.get_all_exons()
    #print 'exon_id', 'exon_stable_id'
    #for e in exons:
    #    print e.id, ' ', e.get_stable_id()
    gene = transcript.get_gene()
    print gene.id, gene.get_stable_id()
    translation = transcript.get_translation()
    print translation.id, translation.get_stable_id()

    geneAdaptor = coreDBAdaptor.get_adaptor('gene')
    gene = geneAdaptor[8946]
    transcripts = gene.get_transcripts()
    for t in transcripts:
        print t.id, t.get_stable_id()

    translationAdaptor = coreDBAdaptor.get_adaptor('translation')
    translation = translationAdaptor[15121]
    transcript = translation.get_transcript()
    print transcript.id, transcript.get_stable_id()
        
    
    
    
    print '\ntest results for the Seqregion class:'
    seq_region = Seqregion(143909)
    print '\nseq_region.getAttributes():'
    seq_region.getAttributes()
    print '\nseq_region.getCoordinateSystem():', seq_region.getCoordinateSystem() # 4
    print '\nseq_region.getName():', seq_region.getName() # AADC01095577.1.1.41877
    print '\nseq_region.getLength():', seq_region.getLength() # 41877
    """
    print '\n\ntest results for the Exon class:'
    exon = driver.get_Adaptor('exon').get_by_dbID(95160)
    #exon = Exon(73777)
    
    #exon = Exon(95172)
    #exon = Exon(95160)
    #print exon.rowobj.seq_region_id # 226034
    #print exon.rowobj.start # 444865
    _sliceable_stableID_tester(exon, 'exon_stable_id')
    """
    _sliceable_tester(exon)
    print '\nmethods unique to the Exon class:'
    print '\nexon.getPhase():', exon.getPhase()
    print '\nexon.getEndPhase():', exon.getEndPhase()
    
    print "\nexon.getSequence('exon')"
    exon_sequence = exon.getSequence('exon')
    print str(exon_sequence)
    print '\nthe length of this exon sequence:', len(exon_sequence)
    
    print '\n\ntest results for the Gene class:'
    #gene = Gene(121)
    gene = Gene(8946)
    _sliceable_stableID_tester(gene, 'gene_stable_id')
    
    _sliceable_tester(gene)
    print '\nmethods unique to the Gene class:'
    print '\ngene.getSequence(\'gene\')'   
    gene_sequence = gene.getSequence('gene')
    print str(gene_sequence)
    print '\nthe length of this gene sequence:', len(gene_sequence)
    print '\ngene.getExons():' 
    transcript_exon_dict = gene.getExons()
    # retrieve and print out all the exons returned by the gene.getExons()    
    if len(transcript_exon_dict) == 0:
        print 'No transcript and therefore no exon identified for this gene.'
    for k, v in transcript_exon_dict.iteritems():
            print '\nall the exons in transcript', k.rowobj.id, ':'
            _getExons_tester(v)        
    print '\ngene.getTranscripts():'
    transcripts = gene.getTranscripts()
    if len(transcripts) == 0:
        print '\nNo transcript identified for this gene.'
    else:
        for index, t in enumerate(transcripts):
            print '\ntranscript ', index, ':'
            t.getAttributes()
            print 'length: ', len(t.getSequence('transcript'))   
    
    # retrieve and print out all the translations returned by the gene.getTranslations()    
    print '\ngene.getTranslations():'
    transcript_translation_dict = gene.getTranslations()
    if len(transcript_translation_dict) == 0:
        print 'No transcript and therefore no translation identified for this gene.'
    for k, v in transcript_translation_dict.iteritems():
            print '\nall the translations in transcript', k.rowobj.id, ':'
            for index, tln in enumerate(v):
                print '\ntranslation', index, ':'
                tln.getAttributes()
    
    print "\ngene.get_created_date('gene_stable_id'):"
    created_date = gene.get_created_date('gene_stable_id')
    print 'created date: ', created_date
    
    
    print '\n\ntest results for the Transcript class:'
    #transcript = Transcript(76)
    transcript = Transcript(15960)
    _sliceable_stableID_tester(transcript, 'transcript_stable_id')
    
    _sliceable_tester(transcript)
    print '\nmethods unique to the Transcript class:'
    print '\ntranscript.getSequence(\'transcript\')'   
    transcript_sequence = transcript.getSequence('transcript')
    print str(transcript_sequence)
    print '\nthe length of this transcript sequence:', len(transcript_sequence)
   
    print '\ntranscript.getExons():' 
    exons = transcript.getExons()
    # retrieve and print out all the exons returned by the transcript.getExons()    
    _getExons_tester(exons)
    
    print '\ntranscript.getTranslations():'
    translations = transcript.getTranslations()
    for index, t in enumerate(translations):
        print '\ntranslation', index, ':'
        t.getAttributes()
        
    # get the gene corresponding to this transcript
    print '\ntranscript.getGene():'
    gene = transcript.getGene()
    print 'The id of its gene:', gene.rowobj.gene_id
    s = gene.getSequence('gene')
    print '\nThe sequence of its gene:', str(s)
    print '\nThe length of its gene:', len(s)
    
    print '\ntest results for the Translation class:'
    translation = Translation(15121)
    print '\ntranslation.getAttributes():'
    translation.getAttributes()
    print '\ntranslation.getTranscriptID:', translation.getTranscriptID()
    print '\ntranslation.get_start_exon_id:', translation.get_start_exon_id()
    print '\ntranslation.get_end_exon_id:', translation.get_end_exon_id()
    print '\ntranslation.getExons():'
    exons = translation.getExons()
    _getExons_tester(exons)
 
    print '\ntest results for the GeneStableID class:'
    gene_stable_id = GeneStableID(8946)
    _stableID_tester(gene_stable_id)

    print '\ntest results for the TranscriptStableID class:'
    transcript_stable_id = TranscriptStableID(15960)
    _stableID_tester(transcript_stable_id)
    
    print '\ntest results for the ExonStableID class:'
    exon_stable_id = ExonStableID(1)
    _stableID_tester(exon_stable_id)

    print '\ntest results for the TranslationStableID class:'
    translation_stable_id = TranslationStableID(1)
    _stableID_tester(translation_stable_id)
    
    print '\ntest results for the Xref class:'
    aXref = Xref(1805202)
    print '\nXref.getAttributes:'
    aXref.getAttributes()
    print '\nXref.get_display_label():', aXref.get_display_label()
    '''
