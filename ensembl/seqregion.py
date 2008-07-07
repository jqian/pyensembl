

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
       # print ('select t1.* from %s t1, %s t2 where t1.cmp_seq_region_id=%%s and t1.asm_seq_region_id=t2.seq_region_id and t2.coord_system_id=%%s and cmp_start-1<=%%s and cmp_end>=%%s' %(self.assembly,self.srdb.seqRegionDB.name)%(srID,self.targetCoord,start,stop))
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
    import MySQLdb
    conn = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous')
    cursor = conn.cursor()
    seq_region = sqlgraph.SQLTable('homo_sapiens_core_47_36i.seq_region',
                                   cursor) 
    dna = sqlgraph.SQLTable('homo_sapiens_core_47_36i.dna', cursor, 
                            itemClass=EnsemblDNA,
                            itemSliceClass=seqdb.SeqDBSlice,
                            attrAlias=dict(seq='sequence'))
    import pygr.Data
    hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
    srdb = SeqRegion(seq_region, {17:hg18, 4:dna}, {17:'chr',4:None})
    chr1 = srdb[226034]
    print len(chr1) # 247249719
    exonSliceDB = sqlgraph.SQLTable('homo_sapiens_core_47_36i.exon',
                                    cursor, itemClass=EnsemblRow)
    #geneSliceDB = sqlgraph.SQLTable('homo_sapiens_core_47_36i.gene',
    #                                cursor, itemClass=EnsemblRow)
    #e = exonSliceDB[73777]
    #print e.start # 444865
    from pygr.seqdb import AnnotationDB
    annoDB = AnnotationDB(exonSliceDB, srdb, sliceAttrDict=
                          dict(id='seq_region_id',stop='seq_region_end',
                               orientation='seq_region_strand'))
    #annoDB = AnnotationDB(geneSliceDB, srdb, sliceAttrDict=
    #                      dict(id='seq_region_id',stop='seq_region_end',
    #                           orientation='seq_region_strand'))
    e = annoDB[73777]
    print e.orientation, e.phase, e.end_phase # 1(always 1!) 1 0
    s = e.sequence
    #g = annoDB[1]
    #print g.orientation, g.is_current, g.analysis_id
    #s = g.sequence
    print len(s)
    print str(s) #'GCAAGCTGTGGACAAGAAGTATGAAGGTCGCTTACAGCATTCTACACAAATTAGGCACAAAGCAGGAACCCATGGTCCGGCCTGGAGATAGG'
    #print str(s.before()[-10:]) # 'GTATTCATAG' last 2 nt is splice-site
    #print str(s.after()[:10]) # 'GTAAGTGCAA' first 2 nt is splice-site
    mapper = EnsemblMapper(annoDB, srdb) # find exons for any sequence slice
    ival = chr1[:100000]
    #ival = chr1[247197525:247197890]
    #ival = chr1[4273:19669]
    print mapper[ival] # find exons in this interval
    #gene_list = mapper[ival]
    #print gene_list[0].orientation
    #s = gene_list[0].sequence
    #print str(s)
    #print len(s)
    print 'neg:', mapper[-ival] # forces annotations to opposite ori...
    #gene_list_neg = mapper[-ival]
    #print gene_list_neg[0].orientation
    #s_neg = gene_list_neg[0].sequence
    #print str(s_neg)
    #print len(s_neg)
    '''
    s = dna[143909] # get this sequence object
    print len(s) # 41877
    print str(s[:10]) # CACCCTGCCC
    amap = AssemblyMapper(srdb, 4, 17) # test the assembly mapper
    contig = srdb[149878] # get a contig
    c = contig[:10] # a short slice of its beginning
    u = amap[c] # map forwards it to hg18
    print str(c) # 'ACCCCTTACC'
    print str(u) # 'accccttacc'
    print c.orientation # 1
    print u.orientation # -1
    ival = chr1[1000000:1000010] # get an hg18 interval
    c = (~amap)[ival] # map back to contig interval
    print repr(ival), repr(c) # chr1[1000000:1000010] 158512[13848:13858]
    print str(ival), str(c) # ACGTGGCTGC ACGTGGCTGC
    '''
    #conn.close()

