# test will be skipped if MySqlDB is unavailable

import unittest, string
from testlib import testutil, logger
from pygr import sqlgraph, seqdb, classutil

class SQLSequence_Test(unittest.TestCase):
    '''Basic SQL sequence class tests
    
    This test setup uses the common (?) method of having the
    SQLSequence objects created by a SQLTable object rather than
    instantiating the SQLSequence objects directly.
    '''
    def setUp(self, cursor=None, dbname='test.sqlsequence_test'):
        createTable = """\
        CREATE TABLE %s
             (primary_id INTEGER PRIMARY KEY %%(AUTO_INCREMENT)s, sequence TEXT)
        """ % dbname
        
        self.db = sqlgraph.SQLTable(dbname, cursor, dropIfExists=True,
            createTable=createTable, attrAlias=dict(seq='sequence'))
        
        self.db.cursor.execute("""\
        INSERT INTO %s (sequence)
              VALUES ('CACCCTGCCCCATCTCCCCAGCCTGGCCCCTCGTGTCTCAGAACCCTCGGGGGGAGGCACAGAAGCCTTCGGGG')
        """ % dbname)

        self.db.cursor.execute("""\
        INSERT INTO %s (sequence)
              VALUES ('GAAAGAAAGAAAGAAAGAAAGAAAGAGAGAGAGAGAGACAGAAG')
        """ % dbname)
        
        class DNASeqRow(seqdb.DNASQLSequence):
            def __len__(self): # just speed optimization
                return self._select('length(sequence)') # SQL SELECT expression
        
        # force the table object to return DNASeqRow objects
        self.db.objclass(DNASeqRow)
        
        self.row1 = self.db[1]
        self.row2 = self.db[2]
        self.EQ = self.assertEqual

    def tearDown(self):
        self.db.cursor.execute('drop table if exists test.sqlsequence_test')

    def test_print(self):
        "Testing identities"
        self.EQ(str(self.row2), 'GAAAGAAAGAAAGAAAGAAAGAAAGAGAGAGAGAGAGACAGAAG')
        self.EQ(repr(self.row2), '2[0:44]')

    def test_len(self):
        "Testing lengths"
        self.EQ(len(self.row2), 44)

    def test_strslice(self):
        "Testing slices"
        self.EQ(self.row2.strslice(3,10), 'AGAAAGA')
    
    def init_subclass_test(self):
        "Testing subclassing"
        self.row2._init_subclass(self.db)

class SQLiteSequence_Test(testutil.SQLite_Mixin, SQLSequence_Test):
    def sqlite_load(self):
        SQLSequence_Test.setUp(self, self.cursor, 'sqlsequence_test')

def get_suite():
    "Returns the testsuite"
    tests = []

    # detect mysql
    if testutil.mysql_enabled():
        tests.append(SQLSequence_Test) 
    else:
        testutil.info('*** skipping SQLSequence_Test')
    if testutil.sqlite_enabled():
        tests.append(SQLiteSequence_Test) 
    else:
        testutil.info('*** skipping SQLSequence_Test')

    return testutil.make_suite(tests)

if __name__ == '__main__':
    suite = get_suite()
    unittest.TextTestRunner(verbosity=2).run(suite)
