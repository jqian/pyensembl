from ensembl.seqregion import *
import MySQLdb



class Driver (object):
    '''Provide access to a database in the core ensembl database system'''
    
    _instance = None

    def __init__(self, host, user, dbname):
        self.conn = MySQLdb.connect(host, user)
        self.db = dbname
        self.tb_adaptor = {'exon': ExonAdaptor, 'gene': GeneAdaptor, 'transcript': TranscriptAdaptor, 'seq_region': SeqregionAdaptor}

    def getAdaptor(self, tbname):        
        adaptor_name = self.tb_adaptor[tbname]
        return adaptor_name(self.db, self.conn.cursor())


def getDriver(host, user, dbname):
        if Driver._instance == None:
            Driver._instance = Driver(host, user, dbname)
        return Driver._instance


class Adaptor(object):
    '''A base class that provides access to a generic table in a core ensembl   
    database'''
   
    def __init__(self, dbname, tbname, RowObj, cursor):
        self.db = dbname
        self.tb = tbname
        self.row = RowObj
        self.tbobj =  sqlgraph.SQLTable(self.db+'.'+self.tb, cursor, 
                                        itemClass=self.row)

    def __getitem__(self, i):
        return self.tbobj[i]


class SeqregionAdaptor(Adaptor):
    '''Provides access to the seq_region table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'seq_region', sqlgraph.TupleO, cursor)

    #def fetch_seqregion_by_something(self, something):

class ExonAdaptor(Adaptor):
    '''Provides access to the exon table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'exon', EnsemblRow, cursor)

    #def fetch_exon_by_something(self, something):


class TranscriptAdaptor(Adaptor):
    'Provides access to the transcript table in an ensembl core database'

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'transcript', EnsemblRow, cursor)

    #def fetch_transcript_by_something(self, something):



class GeneAdaptor(Adaptor):
    '''Provides access to the gene table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'gene', EnsemblRow, cursor)

    #def fetch_genes_by_something(self, something):
    
   
    
if __name__ == '__main__': # example code
    driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
    exontb = driver.getAdaptor('exon')
    exon = exontb[73777]
    print exon.seq_region_id # 226034
    print exon.start # 444865



