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
    '''An interface to a generic record in a table that has a slice of sequence associated to it based on its seq_region_start, seq_region_end and seq_region_strand values'''

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
       
        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        seq_regiontb = driver.getAdaptor('seq_region').tbobj
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        srdb = SeqRegion(seq_regiontb, {17:hg18}, {17:'chr'})
        return srdb 

    def _getAnnotationDB(self, unit_tbname):
        
        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        unit_TB = driver.getAdaptor(unit_tbname).tbobj
        srdb =self._getSeqregionDB()
        
        #for testing purposes
        #print 'seq_region_id =', self.rowobj.seq_region_id # 226034
        #print len(aSeqregion) # 247249719 

        #for testing purposes
        #exon = Exon(73777)   
        #print exon.rowobj.start
        #print exon.rowobj.seq_region_start

        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(unit_TB, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))
        return annoDB
    
    def getSeqregionSeq(self):
        
        srdb = self._getSeqregionDB()
        seq_region_ID = self.getSeqregionID()
        aSeqregion = srdb[seq_region_ID]
        aStart = self.rowobj.start 
        aEnd = self.getSeqregionEnd() - 1
        aSequence = aSeqregion[aStart:aEnd]
        return aSequence
    
    def getSequence(self, unit_tbname):
        'get a sequence object of a sliceable object'

        '''
        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        unit_TB = driver.getAdaptor(unit_tbname).tbobj
        srdb =self._getSeqregionDB()
        
        #for testing purposes
        #print 'seq_region_id =', self.rowobj.seq_region_id # 226034
        #print len(aSeqregion) # 247249719 

        #for testing purposes
        #exon = Exon(73777)   
        #print exon.rowobj.start
        #print exon.rowobj.seq_region_start

        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(unit_TB, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))
        '''
        annoDB = self._getAnnotationDB(unit_tbname)
        
        unitID = self.rowobj.id
        orientation = self.getOrientation()
        unitobj = annoDB[unitID]
        print '+unitobj.id:', unitobj.id
        print '+unitobj.orientation:', unitobj.orientation
        print '+unitobj.start:', unitobj.start
        print '+unitobj.stop:', unitobj.stop
        s = unitobj.sequence
        ssubfirst = s[:10]
        ssublast = s[-10:]
        print '+unitobj.sequence.first10:', str(ssubfirst)
        print '+unitobj.sequence.last10:', str(ssublast)
        if orientation == -1:
            unitobj = -annoDB[unitID]
            
            print '-unitobj.id:', unitobj.id
            print '-unitobj.orientation:', unitobj.orientation
            print '-unitobj.start:', unitobj.start
            print '-unitobj.stop:', unitobj.stop
            s = unitobj.sequence
            ssubfirst = s[:10]
            ssublast = s[-10:]
            print '-unitobj.sequence.first10:', str(ssubfirst)
            print '-unitobj.sequence.last10:', str(ssublast)
        
        s = unitobj.sequence
            
        print 'unitobj.id:', unitobj.id
        print 'unitobj.orientation:', unitobj.orientation
        print 'unitobj.start:', unitobj.start
        print 'unitobj.stop:', unitobj.stop
        ssubfirst = s[:10]
        ssublast = s[-10:]
        print 'unitobj.sequence.first10:', str(ssubfirst)
        print 'unitobj.sequence.last10:', str(ssublast)
        
        print 'length of the sequence: ', len(s)
        print str(s)
        
        return s
        
    def getUnits(self, unit_tbname):
        '''
        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        unit_TB = driver.getAdaptor(unit_tbname).tbobj
        srdb =self._getSeqregionDB(driver)
        
        #for testing purposes
        #print 'seq_region_id =', self.rowobj.seq_region_id # 226034
        #print len(aSeqregion) # 247249719 

        #for testing purposes
        #exon = Exon(73777)   
        #print exon.rowobj.start
        #print exon.rowobj.seq_region_start

        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(unit_TB, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))
        '''
        
        srdb =self._getSeqregionDB()
        annoDB = self._getAnnotationDB(unit_tbname)
        #e = annoDB[73777]
        #print e.phase, e.end_phase
        #s = e.sequence
        #print str(s) # test

        mapper = EnsemblMapper(annoDB, srdb) # find exons for any sequence slice
        ival = self.getSeqregionSeq()
        units = mapper[ival] # find smaller units in this sequence interval
        return units

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

    def getUnits(self):
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

    #def getExons(self):

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

    def getGenes(self):
        'get all the exons of this gene'

        gene_list = self.getUnits('gene')
        print 'genes in interval', self.rowobj.start, '-', self.rowobj.seq_region_end-1, ':'
        print gene_list # find exons in this interval
        genes = []
        for annobj in gene_list:
            g = Gene(annobj.id)
            genes.append(g)

        ''' driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        seq_regiontb = driver.getAdaptor('seq_region').tbobj
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        srdb = SeqRegion(seq_regiontb, {17:hg18}, {17:'chr'})
        aSeqregion = srdb[self.rowobj.seq_region_id]
        
        #for testing purposes
        #print 'seq_region_id =', self.rowobj.seq_region_id # 226034
        #print len(aSeqregion) # 247249719  

        exontb = driver.getAdaptor('exon').tbobj
                    
        #for testing purposes
        #exon = Exon(73777)   
        #print exon.rowobj.start
        #print exon.rowobj.seq_region_start

        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(exontb, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))
        
        #e = annoDB[73777]
        #print e.phase, e.end_phase
        #s = e.sequence
        #print str(s) # test

        mapper = EnsemblMapper(annoDB, srdb) # find exons for any sequence slice
        ival = aSeqregion[self.rowobj.start:(self.rowobj.seq_region_end-1)]
        #ival = aSeqregion[:100000]
        print 'exons in interval', self.rowobj.start, '-', self.rowobj.seq_region_end-1, ':'
        print mapper[ival] # find exons in this interval
        
        exon_list = mapper[ival]
        exons = []
        for annobj in exon_list:
            e = Exon(annobj.id)
            exons.append(e)
        '''  
        for index, g in enumerate(genes):
            print '\ngene ', index, ':'
            g.getAttributes()
            #print 'Phase:', e.getPhase()
            #print 'EndPhase:', e.getEndPhase()
            #print 'Start:', e.rowobj.start

        #e1 = exon_ltist[0]
        #s = e1.id
        #print e1.id
        
        #test
        #e = annoDB[277450]
        #s = e.sequence
        #print e.id
        #print len(s)
        #slice = s[0:365]
        #print str(slice)

    def getExons(self):
        'get all the exons of this gene'

        exon_list = self.getUnits('exon')
        print 'exons in interval', self.rowobj.start, '-', self.rowobj.seq_region_end-1, ':'
        print exon_list # find exons in this interval
        exons = []
        for annobj in exon_list:
            e = Exon(annobj.id)
            exons.append(e)

        ''' driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        seq_regiontb = driver.getAdaptor('seq_region').tbobj
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        srdb = SeqRegion(seq_regiontb, {17:hg18}, {17:'chr'})
        aSeqregion = srdb[self.rowobj.seq_region_id]
        
        #for testing purposes
        #print 'seq_region_id =', self.rowobj.seq_region_id # 226034
        #print len(aSeqregion) # 247249719  

        exontb = driver.getAdaptor('exon').tbobj
                    
        #for testing purposes
        #exon = Exon(73777)   
        #print exon.rowobj.start
        #print exon.rowobj.seq_region_start

        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(exontb, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))
        
        #e = annoDB[73777]
        #print e.phase, e.end_phase
        #s = e.sequence
        #print str(s) # test

        mapper = EnsemblMapper(annoDB, srdb) # find exons for any sequence slice
        ival = aSeqregion[self.rowobj.start:(self.rowobj.seq_region_end-1)]
        #ival = aSeqregion[:100000]
        print 'exons in interval', self.rowobj.start, '-', self.rowobj.seq_region_end-1, ':'
        print mapper[ival] # find exons in this interval
        
        exon_list = mapper[ival]
        exons = []
        for annobj in exon_list:
            e = Exon(annobj.id)
            exons.append(e)
        '''  
        for index, e in enumerate(exons):
            print '\nexon ', index, ':'
            e.getAttributes()
            print 'Phase:', e.getPhase()
            print 'EndPhase:', e.getEndPhase()
            print 'Start:', e.rowobj.start

        #e1 = exon_ltist[0]
        #s = e1.id
        #print e1.id
        
        #test
        #e = annoDB[277450]
        #s = e.sequence
        #print e.id
        #print len(s)
        #slice = s[0:365]
        #print str(slice)

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
    '''A private method to test all the functions defined in the Sliceable class when called by its subclass object.'''

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
    exon = Exon(12)
    #exon = Exon(73777)
    #print exon.rowobj.seq_region_id # 226034
    #print exon.rowobj.start # 444865
    _sliceable_tester(exon)
    print '\nmethods unique to the Exon class:'
    print '\nexon.getPhase():', exon.getPhase()
    print '\nexon.getEndPhase():', exon.getEndPhase()
    exon.getSequence('exon')
   
    print '\n\ntest results for the Gene class:'
    #gene = Gene(24573)
    gene = Gene(1)
    _sliceable_tester(gene)
    print '\nmethods unique to the Gene class:'
    print '\ngene.getExons():'
    gene.getExons()
    gene.getGenes()
    #gene.getSequence('gene')

    print '\n\ntest results for the Transcript class:'
    transcript = Transcript(5)
    _sliceable_tester(transcript)
    #print '\nmethods unique to the Transcript class:'
