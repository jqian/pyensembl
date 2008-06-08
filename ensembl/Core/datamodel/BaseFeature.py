from pygr import sqlgraph
from seqregion import *

class BaseModel(object):
    '''
    '''
    def __init__(self, dbname, tbname, RowObj, cursor):
        self.db = dbname
        self.tb = tbname
        self.row = RowObj
        #self.tbobj =  sqlgraph.SQLTable(self.db+'.'+self.tb, cursor, 
        #                                itemClass=self.row)

        self.tbobj = sqlgraph.SQLTable(self.db+'.'+self.tb, cursor, itemClass=EnsemblRow)
    def __getitem__(self, i):
        return self.tbobj[i]

    def getAttributes(self, i):
        'print out ith row record'
	for k, v in self.tbobj[i]._attrcol.iteritems():
	    print k, ' = ', self.tbobj[i].data[v]
