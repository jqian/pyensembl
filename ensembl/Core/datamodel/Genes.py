from BaseFeature import *
from seqregion import *
from pygr import sqlgraph


class GeneTable(BaseModel):
    'an interface to the gene table in an ensembl core database'

    def __init__(self, dbname, cursor):
        BaseModel.__init__(self, dbname, 'gene', EnsemblRow, cursor)

    def getExons(self, i):
        'get all the exons of the gene record with the primary key i'
       
        from Seqregions import SeqregionTable
        seq_regiondb = SeqregionTable('homo_sapiens_core_47_36i', cursor)
        
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        srdb = SeqRegion(seq_regiondb, {17:hg18}, {17:'chr'})
        aSeqregion = srdb[self.tbobj[i].seq_region_id]
        
        #for testing purposes
        print 'seq_region_id =', self.tbobj[i].seq_region_id # 226034
        print len(aSeqregion) # 247249719        
        
        from Exons import ExonTable
        exondb = ExonTable('homo_sapiens_core_47_36i', cursor)
        exon = exondb[73777]
        
        #for testing purposes
        print exon.start
        print exon.seq_region_start

        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(exondb, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))
        
        e = annoDB[73777]
        print e.phase, e.end_phase # test
        s = e.sequence
        print str(s) # test

        mapper = EnsemblMapper(annoDB, srdb) # find exons for any sequence slice
        ival = aSeqregion[self.tbobj[i].start:self.tbobj[i].end]
        print mapper[ival] # find exons in this interval

    

if __name__ == '__main__': # example code
    import MySQLdb
    conn = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous')
    cursor = conn.cursor()
    
    '''Inclusion of the following block of code invalidates the genedb.getAttrib    utes() method.
    from Exons import ExonTable
    exondb = ExonTable('homo_sapiens_core_47_36i', cursor)
    exondb.getAttributes(73777)
    exon = exondb[73777]
    print exon.start
    print exon.seq_region_start
    '''
    genedb = GeneTable('homo_sapiens_core_47_36i', cursor)
    genedb.getAttributes(2)
    # genedb = sqlgraph.SQLTable('homo_sapiens_core_47_36i.gene',
    #                                cursor, itemClass=EnsemblRow) 
    
    '''Testing:
    gene = genedb.tbobj[2]
    print gene._attrcol
    print gene.data
    '''
    gene = genedb[2]
    print gene.status #KNOWN
    
    '''
    # TESTING Chris'code
    seq_region = sqlgraph.SQLTable('homo_sapiens_core_47_36i.seq_region',
                                   cursor)
    
    exonSliceDB = sqlgraph.SQLTable('homo_sapiens_core_47_36i.exon',
                                    cursor, itemClass=EnsemblRow)
    e = exonSliceDB[73777]
    print e.start
    '''
    #genedb.getExons(24573)
    '''exoncursor = conn.cursor()
    from Exons import ExonTable
    exondb = ExonTable('homo_sapiens_core_47_36i', exoncursor)
    exon = exondb[73777]
    print exon.start
    print exon.seq_region_start
    '''
    
    '''Testing:
    from Seqregions import SeqregionTable
    seq_regiondb = SeqregionTable('homo_sapiens_core_47_36i', cursor)
    seq_region = seq_regiondb[143909]
    print seq_region.name
    '''
        
