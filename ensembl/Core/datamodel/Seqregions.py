from BaseFeature import *
from seqregion import *
from pygr import sqlgraph

class SeqregionTable(BaseModel):
    'an interface to the seq_region table in any ensembl Core database'

    def __init__(self, dbname, cursor):
        BaseModel.__init__(self, dbname, 'seq_region', sqlgraph.TupleO, cursor)
        

    
if __name__ == '__main__': # example code
    import MySQLdb
    conn = MySQLdb.connect(host='ensembldb.ensembl.org', user='anonymous')
    cursor = conn.cursor()
    seqregiondb = SeqregionTable('homo_sapiens_core_47_36i', cursor)
    seq_region = seqregiondb[143909]
    print seq_region.name #AADC01095577.1.1.41877
    seqregiondb.getAttributes(143909)
