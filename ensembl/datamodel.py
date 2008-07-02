from ensembl.adaptor import *



class BaseModel(object):
    '''A generic interface to a row object in a table in a core ensembl database
    '''
    
    def __init__(self, tbname, i):
        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        self.rowobj = driver.getAdaptor(tbname)[i] 

    def getAttributes(self):
        'print out this row record'
	for k, v in self.rowobj._attrcol.iteritems():
	    print k, ' = ', self.rowobj.data[v]



class Seqregion(BaseModel):
    '''An interface to a seq_region record in the seq_region table in any 
    ensembl core database'''

    def __init__(self, i):
        BaseModel.__init__(self, 'seq_region', i)

    def getCoordinateSystem(self):
        return self.rowobj.coord_system_id

    def getName(self):
        return self.rowobj.name

    def getLength(self):
        return self.rowobj.length


class Sliceable(BaseModel):
    '''An interface to a generic record in a table that has a seq_region assigneto it based on the ensembl coordinate system (seq_region_id, seq_region_start, seq_region_end and seq_region_strand)'''

    def __init__(self, tbname, i):
        BaseModel.__init__(self, tbname, i)

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

        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        seq_regiontb = driver.getAdaptor('seq_region').tbobj
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        srdb = SeqRegion(seq_regiontb, {17:hg18}, {17:'chr'})
        return srdb 


    def _getAnnotationDB(self, unit_tbname):
        'a private helper method to create an annotation db object'

        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        unit_TB = driver.getAdaptor(unit_tbname).tbobj
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
        
    def getExons(self):
        'Find exons of a gene or a transcript object.'
        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
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

    

class Exon(Sliceable):
    '''An interface to an exon record in the exon table in an ensembl core 
    database'''

    def __init__(self, i):
        Sliceable.__init__(self, 'exon', i)
    
    def getPhase(self):
        return self.rowobj.phase
    
    def getEndPhase(self):
        return self.rowobj.end_phase

    def getAnalysisID(self):
        return 'None'

    def getXrefID(self):
        return 'None'

    def getBiotype(self):
        return 'None'

    def getStatus(self):
        return 'None'

    def getDescription(self):
        return 'None'

    def isKnown(self):
        return 'Undefined'

    def getExons(self):
        return 'Undefined'

    #def getGene(self):

    #def getTranscripts(self):


class Transcript(Sliceable):
    '''An interface to a transcript record in the transcript table in an ensembl core database'''

    def __init__(self, i):
        Sliceable.__init__(self, 'transcript', i)

    #def getExternalRefs(self):
    #    return getExternalRefs(true);

    #def getExternalRefs(self, includeTranslation):

    #def getCreatedDate(self):

    #def getModifiedDate(self):

    #def getGene(self):

    #def getTranslation(self):

    #def getSequence(self):

    #def getSupportingFeatures(self):

    #def getAnalysis(self):

    #def getGene(self):


class Gene(Sliceable):
    '''An interface to a gene record in the gene table in an ensembl core database'''

    def __init__(self, i):
        Sliceable.__init__(self, 'gene', i)

    def getTranscripts(self):
        'return transcripts if available, otherwise empty list'

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
    
    
if __name__ == '__main__': # example code
    
    print '\ntest results for the Seqregion class:'
    seq_region = Seqregion(143909)
    print '\nseq_region.getAttributes():'
    seq_region.getAttributes()
    print '\nseq_region.getCoordinateSystem():', seq_region.getCoordinateSystem() # 4
    print '\nseq_region.getName():', seq_region.getName() # AADC01095577.1.1.41877
    print '\nseq_region.getLength():', seq_region.getLength() # 41877
    
    print '\n\ntest results for the Exon class:'
    exon = Exon(73777)
    #print exon.rowobj.seq_region_id # 226034
    #print exon.rowobj.start # 444865
    _sliceable_tester(exon)
    print '\nmethods unique to the Exon class:'
    print '\nexon.getPhase():', exon.getPhase()
    print '\nexon.getEndPhase():', exon.getEndPhase()
    print "\nexon.getSequence('exon')"
    exon_sequence = exon.getSequence('exon')
    print str(exon_sequence)
    print '\nthe length of this exon sequence:', len(exon_sequence)
    print '\nexon.getExons():', exon.getExons()

    print '\n\ntest results for the Gene class:'
    gene = Gene(34)
    _sliceable_tester(gene)
    print '\nmethods unique to the Gene class:'
    print '\ngene.getSequence(\'gene\')'   
    gene_sequence = gene.getSequence('gene')
    print str(gene_sequence)
    print '\nthe length of this gene sequence:', len(gene_sequence)
    print '\ngene.getExons():' 
    exons = gene.getExons()
    # retrieve and print out all the exons returned by the gene.getExons()    
    for index, e in enumerate(exons):
            print '\nexon ', index, ':'
            e.getAttributes()
            print 'Phase:', e.getPhase()
            print 'EndPhase:', e.getEndPhase()
            #print 'Sequence:', str(e.getSequence())
            print 'Length of the sequence:', len(e.getSequence('exon'))
    
    print '\n\ntest results for the Transcript class:'
    transcript = Transcript(76)
    _sliceable_tester(transcript)
    print '\nmethods unique to the Transcript class:'
    print '\ntranscript.getSequence(\'transcript\')'   
    transcript_sequence = transcript.getSequence('transcript')
    print str(transcript_sequence)
    print '\nthe length of this transcript sequence:', len(transcript_sequence)
    print '\ntranscript.getExons():' 
    exons = transcript.getExons()
    # retrieve and print out all the exons returned by the transcript.getExons()    
    for index, e in enumerate(exons):
            print '\nexon ', index, ':'
            e.getAttributes()
            print 'Phase:', e.getPhase()
            print 'EndPhase:', e.getEndPhase()
            #print 'Sequence:', str(e.getSequence())
            print 'Length of the sequence:', len(e.getSequence('exon'))
    
