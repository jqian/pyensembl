from ensembl.adaptor import *
from ensembl.seqregion import *
from seqregion import EnsemblRow

def _getDriver():
    '''Obtain an ensembl database connection.  
    Thus, in order to connect to a different database, this is the only place needs to be modified!'
    '''
    
    # Weird...but Python somehow "forgot" the definition and location of getDriver
    from ensembl.adaptor import getDriver
    return getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')

mapperResourceID = 'Bio.Mapping.EnsemblMapper'


class BaseModel(sqlgraph.TupleO):
    '''A generic interface to an item object in a table in the ensembl database
    '''
    
        
    def getAttributes(self):
        'print out this row record'
	#for k, v in self.rowobj._attrcol.iteritems():
        #    print k, ' = ', self.rowobj.data[v]
        for k, v in self._attrcol.iteritems():
	    print k, ' = ', self.data[v]



class Xref(BaseModel):
    '''An interface to a record in the xref table in any ensembl core database'''

    
    #def get_display_label(self):
    #    return self.rowobj.display_label



class Translation(BaseModel):
    '''An interface to an item in the translation table in any ensembl core database'''

    def getExons(self):
        transcript_id = self.rowobj.transcript_id
        start_exon_id = self.get_start_exon_id()
        end_exon_id = self.get_end_exon_id()
        driver = self.driver
        exon_adaptor = driver.getAdaptor('exon')
        exons = exon_adaptor.fetch_exons_by_translation(transcript_id, start_exon_id, end_exon_id)
        return exons
    '''
    def getCreatedDate(self):
        'Get ensembl created date for this translation'

        translation_stableID_adaptor = self.driver.getAdaptor('translation_stable_id')
        created_date = translation_stableID_adaptor[self.translation_id].rowobj.
    '''

class StableID(BaseModel):
    '''An interface to a generic stable_id record in a generic stable_id table in any ensembl core database'''


    #def getStableID(self):
    #    return self.rowobj.stable_id

    #def getVersion(self):
    #    return self.rowobj.version

    #def get_created_date(self):
    #    return self.rowobj.created_date

    #def get_modified_date(self):
    #    return self.rowobj.modified_date


class GeneStableID(StableID):
    '''An interface to a record in the gene_stable_id table in any ensembl core database'''


class TranscriptStableID(StableID):
    '''An interface to a record in the transcript_stable_id table in any ensembl core database'''

class ExonStableID(StableID):
    '''An interface to a record in the exon_stable_id table in any ensembl core database'''

class TranslationStableID(StableID):
    '''An interface to a record in the translation_stable_id table in any ensembl core database'''

class PeptideArchive(BaseModel):
    '''
    An interface to a row record in the peptide_archive table
    >>> driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
    >>> peptide_archive_adaptor = driver.get_Adaptor('peptide_archive')
    >>> peptide_archive = peptide_archive_adaptor.get_by_dbID(100)
    >>> peptide_archive.getAttributes()
    peptide_seq  =  MKKARNDEYENLFNMIVEIPRWTNAKMEIATKEPMNPIKQYVKDGKLRYVANIFPYKGYIWNYGTLPQTWEDPHEKDKSTNCFGDNDPIDVCEIGSKILSCGEVIHVKILGILALIDEGETDWKLIAINANDPEASKFHDIDDVKKFKPGYLEATLNWFRLYKVPDGKPENQFAFNGEFKNKAFALEVIKSTHQCWKALLMKKCNGGAINCTNVQISDSPFRCTQEEARSLVESVSSSPNKESNEEEQVWHFLGK
    peptide_archive_id  =  100
    id  =  100
    md5_checksum  =  A9E4359D28F51F9FF317B378C168BF8D
    '''


"""
class GeneArchive(BaseModel):
    '''
    An interface to a row record in the gene_archive table

    
    '''

    def __init__(self, rowobj):
        BaseModel.__init__(self, rowobj)
"""

class PredictionExon(BaseModel):
    '''
    An interface to a row record in the prediction_exon table

    >>> driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
    >>> predic_exon_adaptor = driver.get_Adaptor('prediction_exon')
    >>> predic_exon = predic_exon_adaptor.get_by_dbID(10)
    >>> predic_exon.getAttributes()
    prediction_exon_id  =  10
    p_value  =  0.039
    seq_region_id  =  149762
    seq_region_end  =  31625
    exon_rank  =  3
    start_phase  =  0
    seq_region_start  =  31474
    score  =  3.84
    prediction_transcript_id  =  3
    id  =  10
    seq_region_strand  =  1
    '''
   
class PredictionTranscript(BaseModel, EnsemblRow):
    '''
    An interface to a prediction_transcript record in any ensembl core database
    
    >>> driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
    >>> predic_transcript_adaptor = driver.get_Adaptor('prediction_transcript')
    >>> predic_transcript = predic_transcript_adaptor.get_by_dbID(150)
    >>> predic_exons = predic_transcript.get_all_Exons()
    >>> for predic_exon in predic_exons:
    ...     print predic_exon.rowobj.prediction_exon_id
    ...
    822
    823
    824
    '''

    def get_all_Exons(self):
        predic_exon_adaptor = self.driver.get_Adaptor('prediction_exon')
        predic_exons = predic_exon_adaptor.get_all_by_ptranscriptID(self.rowobj.prediction_transcript_id)
        return predic_exons


class Seqregion(BaseModel):
    '''An interface to a seq_region record in the seq_region table in any 
    ensembl core database'''


    #def getCoordinateSystem(self):
    #    return self.rowobj.coord_system_id

    #def getName(self):
    #    return self.rowobj.name

    #def getLength(self):
    #    return self.rowobj.length


class Sliceable(BaseModel):
    '''An interface to a generic record in a table that has a seq_region assigneto it based on the ensembl coordinate system (seq_region_id, seq_region_start, seq_region_end and seq_region_strand)'''

    def __init__(self, rowobj):
        BaseModel.__init__(self, rowobj)

    def getSeqregionID(self):
        return self.rowobj.seq_region_id

    def getSeqregionStart(self):
        return self.rowobj.seq_region_start

    def getSeqregionEnd(self):
        return self.rowobj.seq_region_end

    def getOrientation(self):
        return self.rowobj.seq_region_strand

    def is_current(self):
        return self.rowobj.is_current
    

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
    '''    
    def getExons(self):
        'Find exons of a gene or a transcript object.'
 
        # obtain a seq_region database object
        srdb =self._getSeqregionDB()
        # obtain an annotation database object
        annoDB = self._getAnnotationDB('exon')
        # find exons for any sequence slice
        mapper = EnsemblMapper(annoDB, srdb)
        # Obtain a sequence interval object based on the seq_region defined for a gene or a transcript record.  The seq_region interval is described by seq_region_id, seq_region_start, seq_region_end and seq_region_strand.
        ival = self._getSeqregionSeq()
        # find exons in this sequence interval
        exon_list = mapper[ival] 

        print 'exons in interval', self.rowobj.seq_region_start, '-', self.rowobj.seq_region_end, 'on both strands:'
        if len(exon_list) == 0:
            print 'None'
        else:
            # print out all the exons in this interval on both strands
            print exon_list
            exons = []  
            for annobj in exon_list:
                # return only exons on the same strand as the gene or transcript
                if annobj.orientation == 1:
                    e = Exon(annobj.id)
                    exons.append(e)
        return exons
    '''    
    def getAnalysisID(self):
        return self.rowobj.analysis_id
        
    def getBiotype(self):
        return self.rowobj.biotype

    def getStatus(self):
        return self.rowobj.status

    def getDescription(self):
        return self.rowobj.description

    def getXrefID(self):
        return self.rowobj.display_xref_id

    def isKnown(self):
        if self.rowobj.status == 'KNOWN':
            return True
        else:
            return False

    def get_stable_id(self, stableID_tb_name):
        '''Return the stable_id of this Sliceable, none if it doesn't have one'''

        stableID_adaptor = self.driver.get_Adaptor(stableID_tb_name)
        sliceable_id = self.rowobj.id
        stableID_record = stableID_adaptor[sliceable_id]
        stable_id = stableID_record.stable_id
        return stable_id

    def get_created_date(self, stableID_tb_name):
        'Obtain the ensembl created date for the Sliceable if it has a stable id'
        stableID_adaptor = self.driver.get_Adaptor(stableID_tb_name)
        sliceable_id = self.rowobj.id
        stableID_record = stableID_adaptor[sliceable_id]
        created_date = stableID_record.created_date
        return created_date
    
    def get_modified_date(self, stableID_tb_name):
        'Obtain the ensembl modified date for the Sliceable if it has a stable id'
        stableID_adaptor = self.driver.get_Adaptor(stableID_tb_name)
        sliceable_id = self.rowobj.id
        stableID_record = stableID_adaptor[sliceable_id]
        modified_date = stableID_record.modified_date
        return modified_date

    def getVersion(self, stableID_tb_name):
        'Obtain the ensembl version for the Sliceable if it has a stable id'
        
        stableID_adaptor = self.driver.get_Adaptor(stableID_tb_name)
        sliceable_id = self.rowobj.id
        stableID_record = stableID_adaptor[sliceable_id]
        version = stableID_record.version
        return version

        

    

class Exon(BaseModel, EnsemblRow):
    '''An interface to an exon record in the exon table in an ensembl core 
    database

    >>> driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
    >>> exon = driver.get_Adaptor('exon').get_by_dbID(95160)
    >>> _sliceable_stableID_tester(exon, 'exon_stable_id')
    <BLANKLINE>
    get_stable_id(stableID_tb_name): ENSE00001493538
    <BLANKLINE>
    get_created_date(stableID_tb_name): 2006-03-10 00:00:00
    <BLANKLINE>
    get_modified_date(stableID_tb_name): 2006-03-10 00:00:00
    <BLANKLINE>
    getVersion(stableID_tb_name): 1
    '''
    #def getGene(self):

    #def getTranscripts(self):


class Transcript(BaseModel, EnsemblRow):
    '''An interface to a transcript record in the transcript table in an ensembl core database'''

    def get_all_exons(self):
        'obtain all the exons of this transcript'

        # Get the transcriptExons mapper from pygr.Data
        resource_id = mapperResourceID + '.transcriptExons'
        #print self.resource_id
        transcriptExons = adaptor._get_resource(resource_id)
        if transcriptExons == None:
            # Create a pygr.Data transcriptExons mapper
            
            transcriptExons = adaptor._create_pygrDataMapper(sourceTBName, targetTBName)
            adaptorClass = CoreDBAdaptor.TBAdaptorClass[tbname]
            rowClass = CoreDBAdaptor.RowClass[tbname]
            #self.tbobj = sqlgraph.SQLTable(self.db+'.'+self.tb, serverInfo=conn)
            #self.tbobj = sqlgraph.SQLTable(self.db+'.'+self.tb, serverInfo=conn,itemClass=self.row)
            tbAdaptor = adaptorClass(dbname, itemClass=rowClass, serverInfo=self.conn)
           
            # Save self.tbobj to pygr.Data
            #_save_resource(resource_id, tbAdaptor)
        #self.cursor = self.tbobj.cursor
        return tbAdaptor
        #return adaptor_name(self.db_species + '_' + self.db_type + '_' + self.db_version, self.conn)transcript_id = self.rowobj.id
        driver = self.driver
        exon_adaptor = driver.get_Adaptor('exon')
        exons = exon_adaptor.fetch_exons_by_transcriptID(transcript_id)
        return exons

    #def getExternalRefs(self):
    #    return getExternalRefs(true);

    #def getExternalRefs(self, includeTranslation):

    #def getCreatedDate(self):

    #def getModifiedDate(self):

    def getGene(self):
        'obtain its gene'

        gene_id = self.getGeneID()
        gene = Gene(gene_id)
        return gene

    def getTranslations(self):
        'obtain its translation'
        
        translation_adaptor = self.driver.get_Adaptor('translation')
        transcript_id = self.rowobj.transcript_id
        translations = translation_adaptor.fetch_translations_by_transcriptID(transcript_id)
        return translations

    #def getSupportingFeatures(self):

    #def getAnalysis(self):


class Gene(Sliceable):
    '''An interface to a gene record in the gene table in an ensembl core database'''

    def __init__(self, rowobj):
        Sliceable.__init__(self, rowobj)
       
    def getTranscripts(self):
        'return transcripts if available, otherwise empty list'

        gene_id = self.rowobj.gene_id
        transcript_adaptor = self.driver.get_Adaptor('transcript')
        transcripts = transcript_adaptor.fetch_transcripts_by_geneID(gene_id)
        return transcripts
    
    def getExons(self):
        'Obtain all the exons for each of its transcripts'

        transcript_exon_dict = {}
        transcripts = self.getTranscripts()
        if len(transcripts) == 0:
            return transcript_exon_dict
        for t in transcripts:
            exons = t.getExons()
            transcript_exon_dict[t] = exons
        return transcript_exon_dict

    def getTranslations(self):
        'Obtain all the translations for each of its transcripts'
        
        transcript_translation_dict = {}
        transcripts = self.getTranscripts()
        if len(transcripts) == 0:
            return transcript_translation_dict
        for t in transcripts:
            translations = t.getTranslations()
            transcript_translation_dict[t] = translations
        return transcript_translation_dict
        
    def getExternalRefs(self):
        '''References to external databases.  They are a collection of all the 
        transcript.externalRefs.
        Returns a list of ExternalRef objects, list is empty if no external 
        references available.
        '''
        return getExternalRefs(self, true)

    def getExternalRefs(self, includeTranscriptsAndTranslations):
        'includeTranscriptsAndTranslations: a boolean flag'

def _sliceable_tester(classobj):
    '''A private helper method to test all the functions defined in the Sliceable class when invoked by a sliceable object.'''

    print '\ngetAttributes():'
    classobj.getAttributes()
    print '\ngetAnalysisID():', classobj.getAnalysisID()
    print '\ngetSeqregionID():', classobj.getSeqregionID()
    print '\ngetSeqregionStart():', classobj.getSeqregionStart()
    print '\ngetSeqregionEnd():', classobj.getSeqregionEnd()
    print '\ngetOrientation():', classobj.getOrientation()
    print '\ngetXrefID():', classobj.getXrefID()
    print '\ngetBiotype():', classobj.getBiotype()
    print '\ngetStatus():', classobj.getStatus()
    print '\ngetDescription():', classobj.getDescription()
    print '\nis_current():', classobj.is_current()
    print '\nisknown():', classobj.isKnown()  
    
def _sliceable_stableID_tester(classobj, stableID_tb_name):
    
    print '\nget_stable_id(stableID_tb_name):', classobj.get_stable_id(stableID_tb_name)
    print '\nget_created_date(stableID_tb_name):', classobj.get_created_date(stableID_tb_name)
    print '\nget_modified_date(stableID_tb_name):', classobj.get_modified_date(stableID_tb_name)
    print '\ngetVersion(stableID_tb_name):', classobj.getVersion(stableID_tb_name)

def _stableID_tester(classobj):
    '''A private helper method to test all the functions defined in the StableID class when invoked by a specialized StableID object.'''

    print '\ngetAttributes():'
    classobj.getAttributes()
    print '\ngetStableID():', classobj.getStableID()
    print '\ngetVersion():', classobj.getVersion()
    print '\nget_created_date():', classobj.get_created_date()
    print '\nget_modified_date():', classobj.get_modified_date()

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
    """
    driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
    predic_transcript_adaptor = driver.get_Adaptor('prediction_transcript')
    predic_transcript = predic_transcript_adaptor.get_by_dbID(150)
    predic_exons = predic_transcript.get_all_Exons()
    for predic_exon in predic_exons:
        print predic_exon.rowobj.prediction_exon_id
    """
    '''
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
