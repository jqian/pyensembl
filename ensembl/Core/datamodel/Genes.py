from BaseFeature import *
from seqregion import *
from pygr import sqlgraph

class GeneTable(BaseModel):
    '''
    '''

    def __init__(self, dbname, cursor):
        BaseModel.__init__(self,dbname,'gene',EnsemblRow, cursor)

    def getExons(self, i):
        'get all the exons of the gene in the ith row of the gene table'
       
        from Seqregions import SeqregionTable
        seq_regiondb = SeqregionTable('homo_sapiens_core_47_36i', cursor)

        from Exons import ExonTable
        exondb = ExonTable('homo_sapiens_core_47_36i', cursor)
        import pygr.Data
        hg18 = pygr.Data.Bio.Seq.Genome.HUMAN.hg18() # human genome
        srdb = SeqRegion(seq_regiondb, {17:hg18}, {17:'chr'})
        aSeqregion = srdb[self.tbobj[i].seq_region_id]

        print len(aSeqregion)
        e = exondb[73777]
        print e.seq_region_start
        print 'seq_region_id =', self.tbobj[i].seq_region_id

        from pygr.seqdb import AnnotationDB
        annoDB = AnnotationDB(exondb, srdb, sliceAttrDict=
                              dict(id='seq_region_id', stop='seq_region_end',
                                   orientation='seq_region_strand'))
        
        e = annoDB[73777]
        print e.phase, e.end_phase 
        s = e.sequence
        print str(s)

        mapper = EnsemblMapper(annoDB, srdb) # find exons for any sequence slice
        ival = aSeqregion[self.tbobj[i].seq_region_start:self.tbobj[i].seq_region_end]
        print mapper[ival] # find exons in this interval

    

if __name__ == '__main__': # example code
    import MySQLdb
    conn = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous')
    cursor = conn.cursor()
    genedb = GeneTable('homo_sapiens_core_47_36i', cursor)
    # genedb = sqlgraph.SQLTable('homo_sapiens_core_47_36i.gene',
    #                                cursor, itemClass=EnsemblRow) 
    #gene = genedb.tbobj[1]
    gene = genedb[2]
    print gene.status #KNOWN
    genedb.getAttributes(2)

    genedb.getExons(24573)
