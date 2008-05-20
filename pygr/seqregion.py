

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
       name = k.id[len(prefix):] # remove the prefix from the true identifier
       sr = self.inverseDB.seqRegionDB.select('where coord_system_id=%s and name=%s',
                                              (coord_system_id,name)).next()
       return sr.id # its seq_region_id


class SeqRegion(dict):
    '''database proxy to multiple seq databases for Ensembl seq_region table,
    which maps seq_region_id to many different databases each with a
    distinct coord_system_id.    
    '''
    def __init__(self, seqRegionDB, coordSystems, prefixDict=None,
                 prefixDefault='chr'):
        '''    seqRegionDB: dict-like interface to seq_region table
    coordSystems: dictionary of the form {coord_system_id:db_name}
    prefixDict: dictionary of the form {coord_system_id:prefix to add
                                                        to sr.name}
    prefixDefault: default prefix to add to sr.name
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
            seqID = self.prefixDict[sr.coord_system_id] + sr.name
        except KeyError:
            seqID = self.prefixDefault + sr.name
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




if __name__ == '__main__': # example code
    import MySQLdb
    conn = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous')
    cursor = conn.cursor()
    seq_region = sqlgraph.SQLTable('homo_sapiens_core_47_36i.seq_region',
                                   cursor) 
    import pygr.Data
    hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
    srdb = SeqRegion(seq_region, {17:hg18})
    chr1 = srdb[226034]
    print len(chr1) # 247249719
    exonSliceDB = sqlgraph.SQLTable('homo_sapiens_core_47_36i.exon',
                                    cursor, itemClass=EnsemblRow)
    e = exonSliceDB[73777]
    print e.start # 444865
    from pygr.seqdb import AnnotationDB
    annoDB = AnnotationDB(exonSliceDB, srdb, sliceAttrDict=
                          dict(id='seq_region_id',stop='seq_region_end',
                               orientation='seq_region_strand'))
    e = annoDB[73777]
    print e.phase, e.end_phase # 1 0
    s = e.sequence
    print str(s) #'GCAAGCTGTGGACAAGAAGTATGAAGGTCGCTTACAGCATTCTACACAAATTAGGCACAAAGCAGGAACCCATGGTCCGGCCTGGAGATAGG'
    print str(s.before()[-10:]) # 'GTATTCATAG' last 2 nt is splice-site
    print str(s.after()[:10]) # 'GTAAGTGCAA' first 2 nt is splice-site
    mapper = EnsemblMapper(annoDB, srdb) # 
    ival = chr1[:100000]
    print mapper[ival] # find exons in this interval
    print 'neg:', mapper[-ival] # forces annotations to opposite ori...

    #conn.close()

