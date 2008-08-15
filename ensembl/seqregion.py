
'''
class SQLTableMinimal:
    "The part of sqlgraph.SQLTable that can be pickled"
    def __init__(self,name,data,description,primary_key):
        self.name=name
        self.data=data
        self.description=description
        self.primary_key=primary_key
'''

class SeqRegionInv(object):
    'inverse mapping seq obj --> Ensembl seq_region_id'
    def __init__(self, inverseDB):
        self.inverseDB = inverseDB
    def __getitem__(self, k):
       coord_system_id = self.inverseDB.dbs[k.db]
       try: # Ensembl doesn't store correct seqID, so find right prefix
           prefix = self.inverseDB.prefixDict[coord_system_id]
       except KeyError:
           prefix = self.inverseDB.prefixDefault
       if prefix is None:
           return k.id
       name = k.id[len(prefix):] # remove the prefix from the identifier
       sr = self.inverseDB.seqRegionDB.select('where coord_system_id=%s and name=%s',
                                              (coord_system_id,name)).next()
       return sr.id # its seq_region_id


class SeqRegion(dict):
    '''database proxy to multiple seq databases for Ensembl seq_region table,
    which maps a seq_region_id to many different databases each with a
    distinct coord_system_id.    
    '''
    def __init__(self, seqRegionDB, coordSystems, prefixDict=None,
                 prefixDefault='chr'):
        '''    seqRegionDB: dict-like interface to seq_region table
    coordSystems: dictionary of the form {coord_system_id:db_name}
    prefixDict: dictionary of the form {coord_system_id:prefix to add
                                                        to sr.name}
    prefixDefault: default prefix to add to sr.name
    NB: a prefix of "" means sr.name will be used unchanged.
    NB: a prefix of None means k.id will be used unchanged, ignoring sr.name.
        '''
        self.seqRegionDB = seqRegionDB
        self.coordSystems = coordSystems
        d = {}
        for c,db in coordSystems.items(): # create a reverse index
            d[db] = c
        self.dbs = d
        if prefixDict is None:
            prefixDict = {}
        self.prefixDict = prefixDict
        self.prefixDefault = prefixDefault
        self.inverseDB = SeqRegionInv(self)
        dict.__init__(self)
    def __getitem__(self,k):
        'return seq object for the specified seq_region_id'
        try: # retrieve from cache
            return dict.__getitem__(self,k)
        except KeyError:
            pass
        sr = self.seqRegionDB[k] # get seq_region info
        try: # get the right sequence database
            genome = self.coordSystems[sr.coord_system_id]
            #print genome
        except KeyError:
            raise KeyError('unknown coordinate system %d' % sr.coord_system_id)
        try: # Ensembl doesn't store correct seqID, so add right prefix
            prefix = self.prefixDict[sr.coord_system_id]
        except KeyError:
            prefix = self.prefixDefault
        if prefix is None:
            seqID = k
        else:
            seqID = prefix + sr.name
            #print seqID
        s = genome[seqID] # get the actual sequence object
        dict.__setitem__(self, k, s) # save in cache
        return s
    def __invert__(self):
        return self.inverseDB
    
        
class SeqRegionStartDescr(object):
    'converts seq_region_start to Python zero-offset coordinate system'
    def __get__(self, obj, objtype):
        return obj.seq_region_start - 1


from pygr import sqlgraph

class EnsemblRow(sqlgraph.TupleO):
    'use this for all Ensembl tables with seq_region_start'
    start = SeqRegionStartDescr()


    

class EnsemblMapper(object):
    def __init__(self, annoDB, seqRegionDB):
        self.annoDB = annoDB
        self.sliceDB = annoDB.sliceDB
        self.seqRegionInv = ~seqRegionDB
    def __getitem__(self, k):
        srID = self.seqRegionInv[k]
        if k.orientation>0:
            start =  k.start
            stop = k.stop
        else:
            start = -k.stop
            stop = -k.start
        l = []
        # need to correct this to handle +1 offsets exactly right!!!
        for s in self.sliceDB.select('where seq_region_id=%s and seq_region_start<%s and seq_region_end>%s',
                                     (srID,stop,start)):
            if s.start < start:
                a_start = start - s.start
            else:
                a_start = 0
            if s.seq_region_end > stop:
                a_stop = stop - s.start
            else:
                a_stop = s.seq_region_end - s.start
            a = self.annoDB[s.id]
            aslice = a[a_start:a_stop]
            if s.seq_region_strand != k.orientation: # on the opposite strand
                aslice = -aslice # force into same orientation as k
            l.append(aslice)
        return l

from pygr import seqdb
class EnsemblDNA(seqdb.DNASQLSequence):
    'interpret row objects as sequence object a la Ensembl dna table'
    def __len__(self): # just speed optimization
        return self._select('length(sequence)') # SQL SELECT expression

class AssemblyMapperInverse(object):
    'inverse of an AssemblyMapper mapping'
    def __init__(self, mapper):
        self.inverseDB = mapper
    def __getitem__(self, k):
        'map to corresponding interval in the source coord system'
        srID = self.inverseDB.seqRegionInv[k]
        start,stop = k._abs_interval
        n = self.inverseDB.cursor.execute(
            'select t1.* from %s t1, %s t2 where t1.asm_seq_region_id=%%s and t1.cmp_seq_region_id=t2.seq_region_id and t2.coord_system_id=%%s and asm_start-1<=%%s and asm_end>=%%s'
            %(self.inverseDB.assembly,
              self.inverseDB.srdb.seqRegionDB.name),
            (srID,self.inverseDB.sourceCoord,start,stop))
        if n!=1:
            raise KeyError('interval out of mapped range or duplicated!')
        t = self.inverseDB.cursor.fetchall()[0]
        s = self.inverseDB.srdb[t[1]][t[4]-1:t[5]] # sequence slice
        if t[6]<0: # reverse orientation mapping
            s = -s
        u = s[start-t[2]+1:stop-t[2]+1]
        if k.orientation<0: # put back oriented with k
            u = -u
        return u
    def __invert__(self):
        return self.inverseDB

class AssemblyMapper(object):
    'test mapping of contig to genomic coord_system'
    def __init__(self, srdb, sourceCoord, targetCoord):
        '''srdb is SeqRegion object for mapping seq_region_id --> seq
        sourceCoord is coord_system_id of origin coords to map from;
        targetCoord is coord_system_id of desired target to map to'''
        self.srdb = srdb
        self.cursor = srdb.seqRegionDB.cursor
        self.assembly = srdb.seqRegionDB.name.split('.')[0] + '.assembly'
        self.sourceCoord = sourceCoord
        self.targetCoord = targetCoord
        self.seqRegionInv = ~srdb
        self.inverseDB = AssemblyMapperInverse(self)
    def __getitem__(self, k):
        'map to corresponding interval in the target coord system'
        srID = k.id #self.seqRegionInv[k]
        start,stop = k._abs_interval
        print 'start:', start, 'stop:', stop
        print ('select t1.* from %s t1, %s t2 where t1.cmp_seq_region_id=%%s and t1.asm_seq_region_id=t2.seq_region_id and t2.coord_system_id=%%s and cmp_start-1<=%%s and cmp_end>=%%s' %(self.assembly,self.srdb.seqRegionDB.name)%(srID,self.targetCoord,start,stop))
        n = self.cursor.execute(
            'select t1.* from %s t1, %s t2 where t1.cmp_seq_region_id=%%s and t1.asm_seq_region_id=t2.seq_region_id and t2.coord_system_id=%%s and cmp_start-1<=%%s and cmp_end>=%%s'
            %(self.assembly,self.srdb.seqRegionDB.name),
            (srID,self.targetCoord,start,stop))
        if n!=1:
            raise KeyError('interval out of mapped range or duplicated!')
        t = self.cursor.fetchall()[0]
        s = self.srdb[t[0]][t[2]-1:t[3]] # sequence slice
        if t[6]<0: # reverse orientation mapping
            s = -s
        u = s[start-t[4]+1:stop-t[4]+1]
        if k.orientation<0: # put back oriented with k
            u = -u
        return u
    def __invert__(self):
        return self.inverseDB



if __name__ == '__main__': # example code

    conn = sqlgraph.DBServerInfo(host='ensembldb.ensembl.org', user='anonymous')
    
    '''
    # test the TranscriptToExon mapper class
    transcriptTB = sqlgraph.SQLTable('homo_sapiens_core_47_36i.transcript', itemClass=EnsemblRow, serverInfo=conn)
    exonTB = sqlgraph.SQLTable ('homo_sapiens_core_47_36i.exon', itemClass=EnsemblRow, serverInfo=conn)
    
    transcriptToExons = TranscriptToExon(transcriptTB, exonTB)
    transcript = transcriptTB[1]
    exons = transcriptToExons[transcript]
    print 'Transcript', transcript.id, 'has', len(exons), 'exons:'
    print 'id', 'start'
    for e in exons:
        print e.id, e.start
    exon = exonTB[1]
    transcripts = (~transcriptToExons)[exon]
    print 'Exon', exon.id, 'belongs to', len(transcripts), 'transcript:'
    print 'id', 'start'
    for t in transcripts:
        print t.id, t.start
    
    # Save the mapping to pygr.Data
    transcriptTB.__doc__ = 'ensembl transcriptTB (human, core, 47_36i)'
    exonTB.__doc__ = 'ensembl exonTB (human, core, 47_36i)'
    TranscriptToExon.__doc__ = 'ensembl mapper transcriptTB -> exonTB (human, core, 47_36i)'

    pygr.Data.Bio.MySQL.EnsemblTB.HUMAN47_36i.transcriptTB = transcriptTB
    pygr.Data.Bio.MySQL.EnsemblTB.HUMAN47_36i.exonTB = exonTB
    pygr.Data.Bio.Mapping.EnsemblMapper.HUMAN47_36i.transcriptExons = transcriptToExons

    # Create and save the transcript -> exons schema relation to pygr.Data
    pygr.Data.schema.Bio.Mapping.EnsemblMapper.HUMAN47_36i.transcriptExons = pygr.Data.ManyToManyRelation(transcriptTB,exonTB,bindAttrs=('exons','transcripts'))
    pygr.Data.save() # SAVE ALL PENDING DATA AND SCHEMA TO RESOURCE DATABASE

    
    example client code:
    geneToExons = pygr.Data.Bio.Mapping.EnsemblMapper.HUMAN47_36i.transcriptExons()
    # Get the set of exons for this transcript
    # Note: Since we bound an exons attribute to each item of genes, we can get the set of exons for a given transcript as easily as transcript.exons.
    for exon in transcript.exons:
        print exon.id, exon.phase

    # Get the set of transcripts for this exon
    mytranscripts = exon.transcripts # Likewise, this is equivalent to mytranscripts = (~transcriptToExons)[exon]
    '''
    
    seq_region = sqlgraph.SQLTable('homo_sapiens_core_47_36i.seq_region',
                                   serverInfo=conn)
    dna = sqlgraph.SQLTable('homo_sapiens_core_47_36i.dna', serverInfo=conn, 
                            itemClass=EnsemblDNA,
                            itemSliceClass=seqdb.SeqDBSlice,
                            attrAlias=dict(seq='sequence'))
    import pygr.Data
    hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
    '''
    seq_region.__doc__ = 'ensembl_seq_region(human, core, 47_36i)'
    pygr.Data.Bio.MySQL.Ensembl.HUMAN.seq_region47_36i = seq_region
    #dna.__doc__ = 'ensembl_dna(human, core, 47_36i)'
    #pygr.Data.Bio.MySQL.Ensembl.HUMAN.dna47_36i = dna
    
    pygr.Data.save()
    '''
    srdb = SeqRegion(seq_region, {17:hg18, 4:dna}, {17:'chr',4:None})
    chr1 = srdb[226034]
    #print len(chr1) # 247249719
    #chr18 = srdb[226035]
    #print len(chr18) # 76117153 
    exonSliceDB = sqlgraph.SQLTable('homo_sapiens_core_47_36i.exon',
                                   serverInfo=conn, itemClass=EnsemblRow)
    #exonSliceDB.__doc__ = 'ensembl_exon(human, core, 47_36i)'
    #pygr.Data.Bio.MySQL.Ensembl.HUMAN.exonSliceDB47_36i = exonSliceDB
    #pygr.Data.save()
    ##e = exonSliceDB[73777]
    ##print e.start # 444865
    #ptSliceDB = sqlgraph.SQLTable('homo_sapiens_core_47_36i.prediction_transcript', cursor, itemClass = EnsemblRow)
    #ptSliceDB.__doc__ = 'Ensembl_prediction_transcripts'
    #res = SQLTableMinimal(ptSliceDB.name,ptSliceDB.data,ptSliceDB.description,ptSliceDB.primary_key)
    #res.__doc__ = 'Ensembl_prediction_transcripts'
    #print "add resource " + str(res)
    #from StringIO import StringIO
    #import pickle,sys
    #src=StringIO()
    #pickler = pickle.Pickler(src)
    #pickler.root=res
    #pickler.sourceIDs={}
    #pickler.dump(res)
    #pygr.Data.Bio.EnsemblPredictionTranscripts = res
    #pygr.Data.save()
    
    #geneSliceDB = sqlgraph.SQLTable('homo_sapiens_core_47_36i.gene',
    #                                cursor, itemClass=EnsemblRow)
    #e = exonSliceDB[73777]
    #print e.start # 444865
    from pygr.seqdb import AnnotationDB
    annoDB = AnnotationDB(exonSliceDB, srdb, sliceAttrDict=
                          dict(id='seq_region_id',stop='seq_region_end',
                               orientation='seq_region_strand'))
    #annoDB = AnnotationDB(ptSliceDB, srdb, sliceAttrDict=
    #                      dict(id='seq_region_id',stop='seq_region_end',
    #                           orientation='seq_region_strand'))
    #annoDB.__doc__ = 'ensembl prediction_transcripts annotation'
    #pygr.Data.Bio.annotations = annoDB
    
    #pygr.Data.save()
    # test
    #for a in annoDB.itervalues():
    #    print 'db, id:', a.db, a.id
    
    #annoDB = AnnotationDB(geneSliceDB, srdb, sliceAttrDict=
    #                      dict(id='seq_region_id',stop='seq_region_end',
    #                           orientation='seq_region_strand'))
    #e = annoDB[73777]
    #print e.phase, e.end_phase # 1 0
    #pt = annoDB[36948]
    #pt_slice = pt.sequence[:100]
    '''
    print 'length of pt_slice:', len(pt_slice)
    print 'pt sequence first 100bp:', str(pt_slice)
    print 'length:', len(pt) # 17784
    print 'db, id, annotationType:', pt.db, pt.id, pt.annotationType
    print 'seq_region_id (wrong):', pt.id
    print 'sequence id:', pt.sequence.id
    print 'sequence stop:', pt.sequence.stop
    print 'sequence start:', pt.sequence.start
    print 'sequence orientation:', pt.sequence.orientation
    print 'seq_region_id (right):', pt.seq_region_id
    print 'prediction_transcript_id:', pt.prediction_transcript_id
    print 'stop (wrong):', pt.stop # 17784
    print 'stop (right):', pt.seq_region_end
    print 'start:', pt.start
    print 'orientation:', pt.orientation
    #print pt.orientation, pt.start, pt.seq_region_start, pt.stop, pt.id, pt.prediction_transcript_id, pt.display_label, pt.seq_region_id
    #s = pt.sequence
    
    s = e.sequence
    #g = annoDB[1]
    #print g.orientation, g.is_current, g.analysis_id
    #s = g.sequence
    #print len(s)
    print str(s) #'GCAAGCTGTGGACAAGAAGTATGAAGGTCGCTTACAGCATTCTACACAAATTAGGCACAAAGCAGGAACCCATGGTCCGGCCTGGAGATAGG'
    print str(s.before()[-10:]) # 'GTATTCATAG' last 2 nt is splice-site
    print str(s.after()[:10]) # 'GTAAGTGCAA' first 2 nt is splice-site
    '''
    mapper = EnsemblMapper(annoDB, srdb) # find exons for any sequence slice
    #mapper.__doc__ = "mapping of an EnsemblPredictionTranscript to a contig sequence interval in human genome 18" 
    #pygr.Data.getResource.addResource('Bio.Ensembl.mapper', mapper)
    #pygr.Data.here.Bio.mapper = mapper
    #import os
    #os.environ['PYGRDATAPATH'] = '.'
    #pygr.Data.here.getResource.addResource('Bio.Ensembl.mapper', mapper)
    #pygr.Data.here.addResource('Bio.Ensembl.mapper', mapper)
    
    #need to save the annoDB and srdb
    #annoDB.__doc__ = "annotation db for ensembl prediction_transcripts"
    #srdb.__doc__ = "seq_region db for ensembl contig and hg18 sequences"
    #pygr.Data.here.Bio.annoDB = annoDB
    #pygr.Data.Bio.annoDB = annoDB
    #pygr.Data.here.srdb = srdb
    #pygr.Data.here.Bio.srdb = srdb
    #pygr.Data.here.mapper = mapper
    #pygr.Data.here.Bio.mapper = mapper
    #pygr.Data.here.schema.mapper = pygr.Data.OneToManyRelation(srdb, annoDB, bindAttrs=('prediction_transcripts', 'genomic_region'))
    #pygr.Data.save()
    
    '''
    pygr.Data.here.annoDB = annoDB
    pygr.Data.here.srdb = srdb
    pygr.Data.here.mapper = mapper
    pygr.Data.here.schema.mapper = pygr.Data.OneToMany(srdb, annoDB, bindAttrs=('prediction_transcripts', 'genomic_region')
    pygr.Data.save()
    
    #pygr.Data.save(layer='here')
    #pygr.Data.here.save()
    
    #ival = chr18[31140100-31141800]
    ival = chr18[31129340:31190706]

    amap = AssemblyMapper(srdb, 4, 17) # test the assembly mapper
    cival = (~amap)[ival] # map back to contig interval
    contig = srdb[149878][pt.sequence.start:pt.sequence.stop]
    annot36948 = amap[contig]
    annot36948_2 = amap[pt.sequence]
    print 'genomic level interval for annot36948:', repr(annot36948)
    print 'genomic level interval for annot36948_2:', repr(annot36948_2)
    print 'contig level interval for annot36948:', repr(pt.sequence)
    print 'genomic level interval:', repr(ival)
    print 'contig level interval:', repr(cival)
    #print repr(ival), repr(cival)
    '''
    ival = chr1[:100000]
    #ival = chr1[247197525:247197890]
    #ival = chr1[4273:19669]
    print mapper[ival] # find exons in this interval
    #print mapper[cival] # find prediction_transcripts in this interval
    #for a in mapper[cival]:
    #    print 'length:', len(a)
    #gene_list = mapper[ival]
    #print gene_list[0].orientation
    #s = gene_list[0].sequence
    #print str(s)
    #print len(s)
    
    print 'neg:', mapper[-ival] # forces annotations to opposite ori...
    #print 'neg:', mapper[-cival] # forces annotations to opposite ori...
    #gene_list_neg = mapper[-ival]
    #print gene_list_neg[0].orientation
    #s_neg = gene_list_neg[0].sequence
    #print str(s_neg)
    #print len(s_neg)
    '''
    s = dna[143909] # get this sequence object
    print len(s) # 41877
    #print str(s[:10]) # CACCCTGCCC
    '''
    amap = AssemblyMapper(srdb, 4, 17) # test the assembly mapper
    contig = srdb[149878] # get a contig
    contig = contig[48848:66632]
    print repr(contig)
    #g = amap[contig]
    #print repr(g)
    '''
    c = contig[:10] # a short slice of its beginning
    u = amap[c] # map forwards it to hg18
    print str(c) # 'ACCCCTTACC'
    print str(u) # 'accccttacc'
    print c.orientation # 1
    print u.orientation # -1
    ival = chr1[1000000:1000010] # get a hg18 interval
    c = (~amap)[ival] # map back to contig interval
    print repr(ival), repr(c) # chr1[1000000:1000010] 158512[13848:13858]
    print str(ival), str(c) # ACGTGGCTGC ACGTGGCTGC
    #conn.close()
    '''   
