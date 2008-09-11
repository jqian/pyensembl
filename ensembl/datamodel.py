import ensembl.adaptor
from ensembl.seqregion import *
from seqregion import EnsemblRow


class BaseModel(sqlgraph.TupleO):
    '''A generic interface to an item object in a table in the ensembl database
    '''
    
    ''' 
    #It's not working any more with the new pygr version!!! How to print out all the colunms of a table now?
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

    def getDB(self):
        return self.db

class MetaCoord(BaseModel):
    '''An interface to a row in the meta_coord table'''

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
    
    >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> pepAdaptor = coreDBAdaptor.get_adaptor('peptide_archive')
    >>> peptide = pepAdaptor[100]
    >>> print peptide.peptide_seq
    MKKARNDEYENLFNMIVEIPRWTNAKMEIATKEPMNPIKQYVKDGKLRYVANIFPYKGYIWNYGTLPQTWEDPHEKDKSTNCFGDNDPIDVCEIGSKILSCGEVIHVKILGILALIDEGETDWKLIAINANDPEASKFHDIDDVKKFKPGYLEATLNWFRLYKVPDGKPENQFAFNGEFKNKAFALEVIKSTHQCWKALLMKKCNGGAINCTNVQISDSPFRCTQEEARSLVESVSSSPNKESNEEEQVWHFLGK
    '''

class Feature(BaseModel, EnsemblRow):
    '''
    A generic interface to a record in an ensembl feature table.  An ensembl feature table contains columns: seq_region_id, seq_region_start, seq_region_end and seq_region_strand
    '''

    def get_sequence(self, flankingSeq=None):
        '''Obtain a sequence object of the given feature, including its flanking region on both sides if required.  Note: if the feature locates on the reverse strand, the sequence returned will be from the reverse strand.

        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> exonAdaptor = coreDBAdaptor.get_adaptor('exon')
        >>> exon = exonAdaptor[73777]
        >>> exonSeq = exon.get_sequence()
        >>> print str(exonSeq), exonSeq.orientation, repr(exonSeq), len(exonSeq)
        GCAAGCTGTGGACAAGAAGTATGAAGGTCGCTTACAGCATTCTACACAAATTAGGCACAAAGCAGGAACCCATGGTCCGGCCTGGAGATAGG -1 -chr10[444865:444957] 92
        >>> exonSlice = exon.get_sequence(10)
        >>> print str(exonSlice), exonSlice.orientation, repr(exonSlice), len(exonSlice)
        GTATTCATAGGCAAGCTGTGGACAAGAAGTATGAAGGTCGCTTACAGCATTCTACACAAATTAGGCACAAAGCAGGAACCCATGGTCCGGCCTGGAGATAGGGTAAGTGCAA -1 -chr10[444855:444967] 112
        '''
  
        tbName = self.db.name.split('.')[1]
        from ensembl.adaptor import _get_db_parameters
        dbParams = _get_db_parameters(self.db.name)
        from ensembl.adaptor import _get_DB_adaptor
        dbAdaptor = _get_DB_adaptor(dbParams[0], dbParams[1], dbParams[2])
        myAnnodb = dbAdaptor._get_annotationDB(tbName, self.db)
        annobj = myAnnodb[self.id]
        mySeq = annobj.sequence
        if flankingSeq is not None:
            mySeq = mySeq.before()[-flankingSeq:] + mySeq + mySeq.after()[:flankingSeq]
        return mySeq
    
    def transform(self, newCoordSystemName):
        from ensembl.adaptor import _get_db_parameters
        dbParams = _get_db_parameters(self.db.name)
        from ensembl.adaptor import _get_DB_adaptor
        dbAdaptor = _get_DB_adaptor(dbParams[0], dbParams[1], dbParams[2])
        amap = dbAdaptor._get_assemblyMapper()
        seq = self.get_sequence()
        if newCoordSystemName == 'chromosome':
            seq = amap[seq]
        if newCoordSystemName == 'contig':
            seq = (~amap)[seq]
        return seq
            
        
class PredictionExon(Feature):
    '''
    An interface to a row record in the prediction_exon table

    >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> pexonAdaptor = coreDBAdaptor.get_adaptor('prediction_exon')
    >>> pExon = pexonAdaptor[10]
    >>> print pExon.start
    31473
    >>> print pExon.seq_region_start
    31474
    >>> print pExon.prediction_transcript_id
    3
    '''

    def get_prediction_transcript(self):
        '''Obtain the prediction transcript the given prediction exon belongs to

        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> pexonAdaptor = coreDBAdaptor.get_adaptor('prediction_exon')
        >>> pExon = pexonAdaptor[10]
        >>> ptranscript = pExon.get_prediction_transcript()
        >>> print ptranscript.id
        3
        '''
        
        ptranscriptPExons = _get_featureMapper('prediction_transcript', 'prediction_exon', self.db.name)
        ptranscript = (~ptranscriptPExons)[self]
        return ptranscript

    
   
class PredictionTranscript(Feature):
    '''
    An interface to a prediction_transcript record in any ensembl core database
    
    >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> ptranscriptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
    >>> ptranscript = ptranscriptAdaptor[150]
    >>> print ptranscript.start
    18007
    >>> print ptranscript.seq_region_start
    18008
    '''

    def get_prediction_exons(self):
        '''Obtain all the prediction exons for this prediction transcript

        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> ptranscriptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
        >>> ptranscript = ptranscriptAdaptor[150]
        >>> pexons = ptranscript.get_prediction_exons()
        >>> for pe in pexons:
        ...     print pe.id
        ...
        822
        823
        824
        '''

        ptranscriptPExons = _get_featureMapper('prediction_transcript', 'prediction_exon', self.db.name)
        pexons = ptranscriptPExons[self]
        return pexons


class Seqregion(BaseModel):
    '''An interface to a seq_region record in the seq_region table in any 
    ensembl core database'''

class StableObj(object):
    '''An interface to a generic rowObj in a table that has an ensembl stable_idassigned to it.  Subclasses of this class are Gene, Transcript, Translation and Exon'''

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

class Exon(StableObj, Feature):
    '''An interface to an exon record in the exon table in an ensembl core 
    database

    >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
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
    '''

    def get_all_transcripts(self):
        '''Obtain all the transcript annotations this exon belongs to

        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> exonAdaptor = coreDBAdaptor.get_adaptor('exon')
        >>> exon = exonAdaptor[28696]
        >>> transcriptAnnots = exon.get_all_transcripts()
        >>> for t in transcriptAnnots:
        ...   print t.id, t.seq_region_start, len(t.sequence)
        5413 2032725 639645
        5429 2032990 637731
        5445 2032990 637731
        5457 2032990 637731
        5460 2032990 637731
        5468 2032990 637731
        5476 2032990 637731
        5491 2032990 638597
        5510 2032990 637731
        5524 2032990 637731
        5534 2032990 637731
        5547 2032990 638597
        5561 2032990 637731
        5575 2032990 637731
        5581 2032990 638597
        5592 2032990 637731
        5607 2032990 637731
        5616 2032990 637731
        5623 2032990 637731
        5630 2094650 582727
        5637 2428398 248979
        5925 2428398 242323
        5930 2483863 186858
        5945 2428398 248979
        '''
        
        exonTranscriptGraphData = _get_featureGraph('exon', 'transcript', self.db.name)
        exonTranscriptGraph = exonTranscriptGraphData[0]
        exonAnnoDB = exonTranscriptGraphData[1]
        exonAnnot = exonAnnoDB[self.id]
        transcriptAnnots = exonTranscriptGraph[exonAnnot]

        return transcriptAnnots

def _get_featureMapper(sourceTBName, targetTBName, name):
    
    from ensembl.adaptor import _get_db_parameters
    dbParams = _get_db_parameters(name)
    dbSpecies = dbParams[0]
    dbType = dbParams[1]
    dbVersion = dbParams[2]
    from ensembl.adaptor import _get_DB_adaptor
    dbAdaptor = _get_DB_adaptor(dbSpecies, dbType, dbVersion)
    mapper = dbAdaptor._fetch_featureMapper(sourceTBName, targetTBName)
    return mapper

def _get_featureGraph(sourceTBName, targetTBName, name):

    from ensembl.adaptor import _get_db_parameters
    dbParams = _get_db_parameters(name)
    dbSpecies = dbParams[0]
    dbType = dbParams[1]
    dbVersion = dbParams[2]
    from ensembl.adaptor import _get_DB_adaptor
    dbAdaptor = _get_DB_adaptor(dbSpecies, dbType, dbVersion)
    graphData = dbAdaptor._fetch_featureGraph(sourceTBName, targetTBName)
    return graphData
    

class Transcript(StableObj, Feature):
    '''An interface to a transcript record in the transcript table in an ensembl core database

    >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
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
    '''
    

    def get_all_exons(self):
        '''obtain all the exon annotations of this transcript

        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
        >>> transcript = transcriptAdaptor[15960]
        >>> exonAnnots = transcript.get_all_exons() 
        >>> for e in exonAnnots:
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

        >>> transcript = transcriptAdaptor[1]
        >>> exonAnnots = transcript.get_all_exons()
        >>> for e in exonAnnots:
        ...   print e.id, len(e.sequence), e.sequence
        1 273 GGAAAGCGGGTCAAGGCGTAGGGCTGGAGGGCAGGGGCGGGCCCTGGGCGTGGGCTGGGGGTCCTGCCCCGGGGCGCACCCCGGGCGAGGGCTGCCCGGAGGAGCCGAGGTTGGCGGACAGCTTGGCCCTGAGCTTGAGGGGAAGGCAGCGATGGGACAAAGGACGGAGGTCTAGGAAGAGGGTCTGCAGAGCAGAAAGCACGGGTAGGGGCGGCCTGACGCTCGGAAGACAACGCATGGGAGCCGTGTGCACGTCGGGAGCTCGGAGTGAGC
        2 155 GCACCATGACTCCTGTGAGGATGCAGCACTCCCTGGCAGGTCAGACCTATGCCGTGCCCTTCATCCAGCCAGACCTGCGGCGAGAGGAGGCCGTCCAGCAGATGGCGGATGCCCTGCAGTACCTGCAGAAGGTCTCTGGAGACATCTTCAGCAGG
        3 99 GTAGAGCAGAGCCGGAGCCAGGTGCAGGCCATTGGAGAGAAGGTCTCCTTGGCCCAGGCCAAGATTGAGAAGATCAAGGGCAGCAAGAAGGCCATCAAG
        4 147 GTGTTCTCCAGTGCCAAGTACCCTGCTCCAGGGCGCCTGCAGGAATATGGCTCCATCTTCACGGGCGCCCAGGACCCTGGCCTGCAGAGACGCCCCCGCCACAGGATCCAGAGCAAGCACCGCCCCCTGGACGAGCGGGCCCTGCAG
        5 141 GAGAAGCTGAAGGACTTTCCTGTGTGCGTGAGCACCAAGCCGGAGCCCGAGGACGATGCAGAAGAGGGACTTGGGGGTCTTCCCAGCAACATCAGCTCTGTCAGCTCCTTGCTGCTCTTCAACACCACCGAGAACCTGTAT
        6 132 AAGAAGTATGTCTTCCTGGACCCCCTGGCTGGTGCTGTAACAAAGACCCATGTGATGCTGGGGGCAGAGACAGAGGAGAAGCTGTTTGATGCCCCCTTGTCCATCAGCAAGAGAGAGCAGCTGGAACAGCAG
        7 198 GTCCCAGAGAACTACTTCTATGTGCCAGACCTGGGCCAGGTGCCTGAGATTGATGTTCCATCCTACCTGCCTGACCTGCCCGGCATTGCCAACGACCTCATGTACATTGCCGACCTGGGCCCCGGCATTGCCCCCTCTGCCCCTGGCACCATTCCAGAACTGCCCACCTTCCACACTGAGGTAGCCGAGCCTCTCAAG
        8 18 ACCTACAAGATGGGGTac
        9 139 acaccacccccaccgcccccaccaccacccccaGCTCCTGAGGTGCTGGCCAGTGCACCCCCACTCCCACCCTCAACCGCGGCCCCTGTAGGCCAAGGCGCCAGGCAGGACGACAGCAGCAGCAGCGCGTCTCCTTCAG
        10 44 TCCAGGGAGCTCCCAGGGAAGTGGTTGACCCCTCCGGTGGCTGG
        11 106 ACTCTGCTAGAGTCCATCCGCCAAGCTGGGGGCATCGGCAAGGCCAAGCTGCGCAGCATGAAGGAGCGAAAGCTGGAGAAGCAGCAGCAGAAGGAGCAGGAGCAAG
        12 39 TGAGAGCCACGAGCCAAGGTGGGCACTTGATGTCGGATC
        13 92 TGCTCCATGGGGGGACGGCTCCACCCAGCCTGCGCCACTGTGTTCTTAAGAGGCTTCCAGAGAAAACGGCACACCAATCAATAAAGAACTGA
        '''
       
        exonTranscriptGraphData = _get_featureGraph('exon', 'transcript', self.db.name)
        exonTranscriptGraph = exonTranscriptGraphData[0]
        transcriptAnnoDB = exonTranscriptGraphData[2]
        transcriptAnnot = transcriptAnnoDB[self.id]
        exonAnnots = (~exonTranscriptGraph)[transcriptAnnot]

        return exonAnnots

    def get_spliced_seq(self):
        '''Obtain the spliced sequence of the transcript
        
        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
        >>> transcript = transcriptAdaptor[1]
        >>> splicedSeq = transcript.get_spliced_seq()
        >>> print len(splicedSeq), str(splicedSeq)
        1583 GGAAAGCGGGTCAAGGCGTAGGGCTGGAGGGCAGGGGCGGGCCCTGGGCGTGGGCTGGGGGTCCTGCCCCGGGGCGCACCCCGGGCGAGGGCTGCCCGGAGGAGCCGAGGTTGGCGGACAGCTTGGCCCTGAGCTTGAGGGGAAGGCAGCGATGGGACAAAGGACGGAGGTCTAGGAAGAGGGTCTGCAGAGCAGAAAGCACGGGTAGGGGCGGCCTGACGCTCGGAAGACAACGCATGGGAGCCGTGTGCACGTCGGGAGCTCGGAGTGAGCGCACCATGACTCCTGTGAGGATGCAGCACTCCCTGGCAGGTCAGACCTATGCCGTGCCCTTCATCCAGCCAGACCTGCGGCGAGAGGAGGCCGTCCAGCAGATGGCGGATGCCCTGCAGTACCTGCAGAAGGTCTCTGGAGACATCTTCAGCAGGGTAGAGCAGAGCCGGAGCCAGGTGCAGGCCATTGGAGAGAAGGTCTCCTTGGCCCAGGCCAAGATTGAGAAGATCAAGGGCAGCAAGAAGGCCATCAAGGTGTTCTCCAGTGCCAAGTACCCTGCTCCAGGGCGCCTGCAGGAATATGGCTCCATCTTCACGGGCGCCCAGGACCCTGGCCTGCAGAGACGCCCCCGCCACAGGATCCAGAGCAAGCACCGCCCCCTGGACGAGCGGGCCCTGCAGGAGAAGCTGAAGGACTTTCCTGTGTGCGTGAGCACCAAGCCGGAGCCCGAGGACGATGCAGAAGAGGGACTTGGGGGTCTTCCCAGCAACATCAGCTCTGTCAGCTCCTTGCTGCTCTTCAACACCACCGAGAACCTGTATAAGAAGTATGTCTTCCTGGACCCCCTGGCTGGTGCTGTAACAAAGACCCATGTGATGCTGGGGGCAGAGACAGAGGAGAAGCTGTTTGATGCCCCCTTGTCCATCAGCAAGAGAGAGCAGCTGGAACAGCAGGTCCCAGAGAACTACTTCTATGTGCCAGACCTGGGCCAGGTGCCTGAGATTGATGTTCCATCCTACCTGCCTGACCTGCCCGGCATTGCCAACGACCTCATGTACATTGCCGACCTGGGCCCCGGCATTGCCCCCTCTGCCCCTGGCACCATTCCAGAACTGCCCACCTTCCACACTGAGGTAGCCGAGCCTCTCAAGACCTACAAGATGGGGTacacaccacccccaccgcccccaccaccacccccaGCTCCTGAGGTGCTGGCCAGTGCACCCCCACTCCCACCCTCAACCGCGGCCCCTGTAGGCCAAGGCGCCAGGCAGGACGACAGCAGCAGCAGCGCGTCTCCTTCAGTCCAGGGAGCTCCCAGGGAAGTGGTTGACCCCTCCGGTGGCTGGACTCTGCTAGAGTCCATCCGCCAAGCTGGGGGCATCGGCAAGGCCAAGCTGCGCAGCATGAAGGAGCGAAAGCTGGAGAAGCAGCAGCAGAAGGAGCAGGAGCAAGTGAGAGCCACGAGCCAAGGTGGGCACTTGATGTCGGATCTGCTCCATGGGGGGACGGCTCCACCCAGCCTGCGCCACTGTGTTCTTAAGAGGCTTCCAGAGAAAACGGCACACCAATCAATAAAGAACTGA
        '''
        
        # stitches all the exon sequences of the transcript together
        exonAnnots = self.get_all_exons()
        splicedSeq = ''
        for e in exonAnnots:
            splicedSeq += str(e.sequence)
        return splicedSeq

    def get_five_utr(self): 
        '''Obtain the five-prime untranslated region of the transcript.  Return None, if the transcript is not translateable or the translation starts at the beginning of the start_exon.

        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
        >>> transcript = transcriptAdaptor[1]
        >>> fiveUtr = transcript.get_five_utr()
        >>> if fiveUtr is not None:
        ...    print len(fiveUtr.sequence), str(fiveUtr.sequence)
        ...
        236 GGAAAGCGGGTCAAGGCGTAGGGCTGGAGGGCAGGGGCGGGCCCTGGGCGTGGGCTGGGGGTCCTGCCCCGGGGCGCACCCCGGGCGAGGGCTGCCCGGAGGAGCCGAGGTTGGCGGACAGCTTGGCCCTGAGCTTGAGGGGAAGGCAGCGATGGGACAAAGGACGGAGGTCTAGGAAGAGGGTCTGCAGAGCAGAAAGCACGGGTAGGGGCGGCCTGACGCTCGGAAGACAACGC
        >>> transcript = transcriptAdaptor[15121]
        >>> fiveUtr = transcript.get_five_utr() # rainy test
        No 5'-untranslated region found.  Translation starts at the beginning of the start_exon
        >>> if fiveUtr is not None:
        ...    print repr(fiveUtr.sequence), len(fiveUtr.sequence), str(fiveUtr.sequence)
        
        >>> threeUtr = transcript.get_three_utr()
        >>> if threeUtr is not None:
        ...    print len(threeUtr.sequence), str(threeUtr.sequence)
        ...
        649 CTACCGTCTTTTTGCTAGGACTTAAACTGACTTGAGTGTGGCAAAAAGTTAACAAAAAAGGAGAAAAAATGAACAATCGTTTGTGGTTTCTTGGGAAAACTTTTCATACCAGGTGATACTATTCAAAAACCCCGTTGTCTCCCTGCAAGTGCTGATTTGAAATGCAGAAGCCACAGTaaaaaaaaaaaaaaaaaaaaaaaaaaagaaaaaaaaaTCAAAATGTATAAATATTGGAAATCAAGTTTTTCAGCTGTTTTGTTGGTTGGTTGGTTGGTTTTTGTTTGGTTTTGTTTAAATGGGCAAGAAGTAAATAATGTGGCTGGAATACAAGTTGAACAAACTAGAAGACACAAATCTAACATAGTTTTTATGGACCAAGGAACTTGTATATTGTATAAGCTTTAGTAAAAGGTACATTTTCACCATACCTTTTTTTATATCACGGTATTATAGTACACCTTGTTACCAAATAGGTTGTTCTCTTCCCCACCCACCTTTGAGCTTTTGCTCTAAAATACATTCAGGTTCCAAGCCTGACCATCCTTGTTTAATCTATCATACTCTTCCAGGTTTTTTTTTTTTGGTCTAAGGCTGGAACTTTTTTCTTTTTTTTCAGCTGAAGTCTTATGACTTTTCATGAGTCAAAATT
        >>> transcript = transcriptAdaptor[214] # boundary test
        >>> fiveUtr = transcript.get_five_utr()
        This transcript is not translateable!
        '''

        fiveUtr = None
        exonAnnots = self.get_all_exons()
        translation = self.get_translation()
        if translation is not None:
            startExonID = translation.start_exon_id
            start = translation.seq_start
            for e in exonAnnots: 
                if e.id == startExonID:
                    startExon = e
                    #print repr(startExon.sequence)
            try:
                fiveUtr = startExon[:start-1]
            except IndexError:
                print "No 5'-untranslated region found.  Translation starts at the beginning of the start_exon"
        return fiveUtr

    def get_three_utr(self):
        '''Obtain the three-prime untranslated region of the transcript.  ReturnNone, if the transcript is not translateable or the translation ends in the end of the end_exon.

        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
        >>> transcript = transcriptAdaptor[1]
        >>> threeUtr = transcript.get_three_utr() # rainy test
        No 3'-untranslated region found.  All the sequence of the end_exon is needed for translation.

        >>> if threeUtr is not None:
        ...    print len(threeUtr.sequence), str(threeUtr.sequence)
        
        >>> transcript = transcriptAdaptor[15121]
        >>> threeUtr = transcript.get_three_utr()
        >>> if threeUtr is not None:
        ...    print len(threeUtr.sequence), str(threeUtr.sequence)
        ...
        649 CTACCGTCTTTTTGCTAGGACTTAAACTGACTTGAGTGTGGCAAAAAGTTAACAAAAAAGGAGAAAAAATGAACAATCGTTTGTGGTTTCTTGGGAAAACTTTTCATACCAGGTGATACTATTCAAAAACCCCGTTGTCTCCCTGCAAGTGCTGATTTGAAATGCAGAAGCCACAGTaaaaaaaaaaaaaaaaaaaaaaaaaaagaaaaaaaaaTCAAAATGTATAAATATTGGAAATCAAGTTTTTCAGCTGTTTTGTTGGTTGGTTGGTTGGTTTTTGTTTGGTTTTGTTTAAATGGGCAAGAAGTAAATAATGTGGCTGGAATACAAGTTGAACAAACTAGAAGACACAAATCTAACATAGTTTTTATGGACCAAGGAACTTGTATATTGTATAAGCTTTAGTAAAAGGTACATTTTCACCATACCTTTTTTTATATCACGGTATTATAGTACACCTTGTTACCAAATAGGTTGTTCTCTTCCCCACCCACCTTTGAGCTTTTGCTCTAAAATACATTCAGGTTCCAAGCCTGACCATCCTTGTTTAATCTATCATACTCTTCCAGGTTTTTTTTTTTTGGTCTAAGGCTGGAACTTTTTTCTTTTTTTTCAGCTGAAGTCTTATGACTTTTCATGAGTCAAAATT
        >>> transcript = transcriptAdaptor[214] # boundary test
        >>> threeUtr = transcript.get_three_utr()
        This transcript is not translateable!
        '''
        
        threeUtr = None
        exonAnnots = self.get_all_exons()
        translation = self.get_translation()
        if translation is not None:
            endExonID = translation.end_exon_id
            end = translation.seq_end
            #print end
            for e in exonAnnots:
                if e.id == endExonID:
                    endExon = e
                    #print e.id, repr(endExon.sequence)
            try:
                threeUtr = endExon[end:]
            except IndexError:
                print "No 3'-untranslated region found.  All the sequence of the end_exon is needed for translation."
        return threeUtr
        

    #def getExternalRefs(self):
    #    return getExternalRefs(true);

    #def getExternalRefs(self, includeTranslation):

    def get_gene(self):
        '''obtain its gene
        
        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
        >>> transcript = transcriptAdaptor[15960]
        >>> gene = transcript.get_gene() 
        >>> print gene.id, gene.get_stable_id()
        8946 ENSG00000120645
        '''
        
        geneTranscripts = _get_featureMapper('gene', 'transcript', self.db.name)
        gene = (~geneTranscripts)[self]
        return gene

    def get_translation(self):
        '''Obtain its translation.  Return None if no translation found for this transcript.
        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
        >>> transcript = transcriptAdaptor[15960]
        >>> translation = transcript.get_translation()
        >>> print translation.id, translation.get_stable_id()
        15121 ENSP00000372292

        >>> transcript = transcriptAdaptor[214] # boundary test
        >>> translation = transcript.get_translation()
        This transcript is not translateable!
        '''
        
        transcriptTranslation = _get_featureMapper('transcript', 'translation', self.db.name)
        translation = transcriptTranslation[self]
        return translation

    #def get_translateable_exons(self):
        'obtain all translateable exons'
        '''
        # get all the exon annotations
        exonAnnots = self.get_all_exons()
        num = len(exonAnnots)
        # get its translation object
        translation = self.get_translation()
        startExonID = translation.start_exon_id
        endExonID = translation.end_exon_id
        # get all the exons starting from the start_exon
        myExonAnnots = []
        count = 0
        myCount = 0
        for e in exonAnnots:
            if e.id == startExonID:
                myExonAnnots.append(e)
                myCount = count
            count = count + 1
        print myCount, count
        l = range(myCount, num)
        print l
        for i in l:
            print i
            if i <= myCount:
                exonAnnots.next()
            else:
                exonAnnot = exonAnnots.next()
                myExonAnnots.append(exonAnnot)
        '''
        '''
        print exonAnnots
        exon = exonAnnots[2] # Is there a way to access one item in the exonAnnots list other than the for-loop?
        print repr(exon)
        
        for index in range(myCount, num):
        
            myExonAnnots.append(exonAnnots[index])
        # get all the translateable exons 
        transExonAnnots = []
        count = 0
        for e in myExonAnnots:
            if e.id != endExonID:
                transExonAnnots.append(e)
                count = count + 1
        print count
        endExonAnnot = myExonAnnots[count-1]
        transExonAnnots.append(endExonAnnot)
        return transExonAnnots
        '''
    #def getSupportingFeatures(self):

    #def getAnalysis(self):


class Gene(StableObj, Feature):
    '''An interface to a gene record in the gene table in an ensembl core database
    '''

    def get_transcripts(self):
        '''return transcripts if available, otherwise empty list
        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
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
        
        geneTranscripts = _get_featureMapper('gene', 'transcript', self.db.name)
        transcripts = geneTranscripts[self]
        return transcripts
        
    #def getExternalRefs(self):
        '''References to external databases.  They are a collection of all the 
        transcript.externalRefs.
        Returns a list of ExternalRef objects, list is empty if no external 
        references available.
        '''
        return getExternalRefs(self, true)

    #def getExternalRefs(self, includeTranscriptsAndTranslations):
        'includeTranscriptsAndTranslations: a boolean flag'


class Translation(StableObj, BaseModel):
    '''An interface to an item in the translation table in any ensembl core database
    '''

    def get_transcript(self):
        '''obtain its transcript
        
        >>> serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> translationAdaptor = coreDBAdaptor.get_adaptor('translation')
        >>> translation = translationAdaptor[15121]
        >>> transcript = translation.get_transcript()
        >>> print transcript.id, transcript.get_stable_id()
        15960 ENST00000382841
        '''
        
        transcriptTranslation = _get_featureMapper('transcript', 'translation', self.db.name)
        transcript = (~transcriptTranslation)[self]
        return transcript
    '''
    def getExons(self):
        transcript_id = self.rowobj.transcript_id
        start_exon_id = self.get_start_exon_id()
        end_exon_id = self.get_end_exon_id()
        driver = self.driver
        exon_adaptor = driver.getAdaptor('exon')
        exons = exon_adaptor.fetch_exons_by_translation(transcript_id, start_exon_id, end_exon_id)
        return exons
    '''

def _test():
    import doctest
    doctest.testmod()

if __name__ == '__main__': # example code
    
    #_test()
    '''
    serverRegistry = ensembl.adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
    coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    
    pexonAdaptor = coreDBAdaptor.get_adaptor('prediction_exon')
    pExon = pexonAdaptor[10]
    print pExon.start
    print pExon.seq_region_start
    print pExon.prediction_transcript_id
    ptranscript = pExon.get_prediction_transcript()
    print ptranscript.id
    
    ptranscriptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
    #ptranscript = ptranscriptAdaptor[3]
    ptranscript = ptranscriptAdaptor[36948]
    #print ptranscript.getDB()
    ptSeq = ptranscript.get_sequence()
    print repr(ptSeq)
    print ptSeq.id, ptSeq.start, ptSeq.stop
    gptSeq = ptranscript.transform('chromosome')
    print repr(gptSeq)
    print gptSeq.id, gptSeq.start, gptSeq.stop
    #print ptranscript.start
    #print ptranscript.seq_region_start
    #pexons = ptranscript.get_prediction_exons()
    #for pe in pexons:
    #    print pe.id
    
    exonAdaptor = coreDBAdaptor.get_adaptor('exon')
    exon = exonAdaptor[28696]
    transcriptAnnots = exon.get_all_transcripts()
    for t in transcriptAnnots:
        print t.id, t.seq_region_start, len(t.sequence)
    
    transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
    transcript = transcriptAdaptor[214]
    #transcript = transcriptAdaptor[1]
    translation = transcript.get_translation()
    if translation is not None:
        print translation.id
    
    #transcript = transcriptAdaptor[1]
    #transcript = transcriptAdaptor[15960]
    #transcript = transcriptAdaptor[15121]
    #exonAnnots = transcript.get_all_exons()
    #for e in exonAnnots:
    #    print e.id, len(e.sequence), e.sequence
    #splicedSeq = transcript.get_spliced_seq()
    #print len(splicedSeq), str(splicedSeq)
    
    fiveUtr = transcript.get_five_utr()
    if fiveUtr is not None:
        print repr(fiveUtr.sequence), len(fiveUtr.sequence), str(fiveUtr.sequence)
    
    threeUtr = transcript.get_three_utr()
    if threeUtr is not None:
        print repr(threeUtr.sequence), len(threeUtr.sequence), str(threeUtr.sequence)
    
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
