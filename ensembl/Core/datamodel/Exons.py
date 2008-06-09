from pygr import sqlgraph
from pygr import seqregion
from seqregion import *
from BaseFeature import *



class ExonTable(BaseModel):
    'an interface to the exon table in an ensembl Core database'

    def __init__(self, dbname, cursor):
        BaseModel.__init__(self, dbname, 'exon', EnsemblRow, cursor)
    
if __name__ == '__main__': # example code
    import MySQLdb
    conn = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous')
    cursor = conn.cursor()
    exondb = ExonTable('homo_sapiens_core_47_36i', cursor)
    exon = exondb[73777]
    print exon.seq_region_id # 226034
    print exon.start # 444865

    exondb.getAttributes(73777)
