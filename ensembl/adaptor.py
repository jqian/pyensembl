import MySQLdb
import pygr.Data
from pygr import seqdb
from pygr import sqlgraph
#from ensembl import seqregion
import seqregion
from seqregion import EnsemblDNA
from seqregion import SeqRegion
from ensembl.featuremapping import *
#from ensembl.datamodel import *
import datamodel
from datamodel import MetaCoord
from datamodel import Dna
from datamodel import Exon
from datamodel import ExonAnnot
from datamodel import Gene
from datamodel import GeneAnnot
from datamodel import Transcript
from datamodel import TranscriptAnnot
from datamodel import Seqregion
from datamodel import Translation
from datamodel import GeneStableID
from datamodel import TranscriptStableID
from datamodel import TranslationStableID
from datamodel import ExonStableID
from datamodel import PredictionTranscript
from datamodel import PredictionExon
from datamodel import PeptideArchive
from datamodel import CoordSystem
from datamodel import Xref

def _get_db_parameters(name):
    dbName = name.split('.')[0]
    dbSpecies = dbName.split('_')[0] + '_' + dbName.split('_')[1]
    dbType = dbName.split('_')[2]
    dbVersion = dbName.split('_')[3] + '_' + dbName.split('_')[4]
    dbParameters = [dbSpecies, dbType, dbVersion]
    return dbParameters



class FeatureAdaptor(sqlgraph.SQLTable):
    '''Provide a generic interface to an ensembl feature table.  Such table contains the 'seq_region_id', 'seq_region_start', 'seq_region_end' and 'seq_region_strand' columns.'''

    def fetch_all_by_slice(self, slice):
        'Find all the features in the given slice'

        # find out slice is at genomic level or contig level
        try:
            sliceName = slice.id[0] + slice.id[1] + slice.id[2]
            #sliceName = slice.id
        except TypeError:
            sliceName = 'None'
        #print sliceName
        if sliceName == 'chr':
            level = 17 # genomic slice
        else:
            level = 4 # contig slice
        # get the DB adaptor
        dbParams = _get_db_parameters(self.name)
        dbAdaptor = _get_DB_adaptor(dbParams[0], dbParams[1], dbParams[2])
        tbName = self.name.split('.')[1]
        metaCoordAdaptor = dbAdaptor.get_adaptor('meta_coord')
        t = metaCoordAdaptor.select('where table_name = %s', (tbName))
        coordSystemID = 17
        for row in t:
            if row.coord_system_id == 4:
                coordSystemID = 4
        #print coordSystemID
        if coordSystemID == 17 and level == 4:
            amap = dbAdaptor._get_assemblyMapper()
            slice = amap[slice] # map to genomic interval
            #print repr(slice)
        if coordSystemID == 4 and level == 17:
            # transform the genomic slice into a contig slice
            amap = dbAdaptor._get_assemblyMapper()
            slice = (~amap)[slice] # map to contig interval
            #print repr(slice)
        #print repr(slice)
        seqFeatureMapper = dbAdaptor._get_seqFeatureMapper(tbName, self)
        features = seqFeatureMapper[slice]
        return features
    
class MetaCoordAdaptor(sqlgraph.SQLTable):
    '''Provides access to the meta_coord table in an ensembl core database'''


class TranslationAdaptor(sqlgraph.SQLTable):
    '''Provides access to the translation table in an ensembl core database

    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> translationAdaptor = coreDBAdaptor.get_adaptor('translation')
    >>> translation = translationAdaptor[10]
    >>> print translation.transcript_id
    10
    >>> print translation.start_exon_id
    61
    >>> print translation.end_exon_id
    62
    >>> translations = translationAdaptor.fetch_by_stable_id('ENSP00000372292')
    >>> for tl in translations:
    ...     print tl.id, tl.get_stable_id()
    ...
    15121 ENSP00000372292

    '''
    def fetch_by_stable_id(self, translationStableID):
        
        dbParams = _get_db_parameters(self.name)
        dbAdaptor = _get_DB_adaptor(dbParams[0], dbParams[1], dbParams[2])
        tlsAdaptor = dbAdaptor.get_adaptor('translation_stable_id')
        t = tlsAdaptor.select('where stable_id = %s', (translationStableID))
        translations = []
        for row in t:
            translation = self[row.translation_id]
            translations.append(translation)
        return translations



def get_registry(**kwargs):
    '''Provides a method to generate a connection object to the MySQL server

    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> exonAdaptor = coreDBAdaptor.get_adaptor('exon')
    >>> exon = exonAdaptor[10]
    >>> print exon.phase
    1
    >>> print exon.seq_region_strand
    -1
    '''
   
    if Registry._instance == None:
        Registry._instance = Registry(**kwargs)
    return Registry._instance



def _get_resource(resource_id):
    '''Obtain the required ensembl table adaptor from pygr.Data, return None if unable to find it in PYGRDATAPATH'''

    import pygr.Data
    try:
        return pygr.Data.getResource(resource_id)
    except pygr.Data.PygrDataNotFoundError:
        return None


def _get_DB_adaptor(dbSpecies, dbType, dbVersion):
    
    serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    myDBAdaptor = serverRegistry.get_DBAdaptor(dbSpecies, dbType, dbVersion)
    return myDBAdaptor


def _save_resource(resource_id, resource, suffix=None):
    'Save an ensembl database resource to pygr.Data'

    import pygr.Data
    if suffix is None:
        resource.__doc__ = 'ensembl ' + resource_id.split('.')[3] + ' ' + resource_id.split('.')[4]
    else:
        resource.__doc__ = 'ensembl ' + resource_id.split('.')[3] + ' ' + resource_id.split('.')[4] + ' ' + suffix 
    print 'saving', resource.__doc__, '...'
    pygr.Data.addResource(resource_id, resource)
    # Save all pending data and schema to resource database
    pygr.Data.save()


#def _save_mapper(mapperID, mapper, sourceDB, targetDB):
#    'Save the ensembl mapper to pygr.Data'


#def _save_graph(graphID, graph, sourceDB, targetDB):
#    'Save the ensembl graph to pygr.Data'

class Registry(object):
    'Provide a connection to the ensembl server'

    _instance = None

    def __init__(self, **kwargs):
        self.kwargs = kwargs # connection arguments
        # get an ensembl server connection from pygr.Data, if it's not there, create one and save it
        connID = 'Bio.Server.Ensembl.Ensembldb'
        self.conn = _get_resource(connID)
        if self.conn is None:
            self.conn = sqlgraph.DBServerInfo(**self.kwargs)
            # save the DBServerInfo for connecting to Ensembl's MySQL server into pygr.Data
            import pygr.Data
            self.conn.__doc__ = 'DBServerInfo for connecting to ensembldb.ensembl.org'
            print 'saving', self.conn.__doc__, '...'
            pygr.Data.addResource(connID, self.conn)
            # Save all pending data and schema to resource database
            pygr.Data.save()
        self.db_adaptor = {'core': CoreDBAdaptor}

    def get_DBAdaptor(self, db_species, db_type, db_version):
        'Obtain an adaptor to the given type of database'

        db_adaptor = self.db_adaptor[db_type]
        return db_adaptor(self.conn, db_species, db_version)


class MetaCoordAdaptor(sqlgraph.SQLTable):
    '''Provides access to the meta_coord table in an ensembl core database'''
    

class XrefAdaptor(sqlgraph.SQLTable):
    '''Provides access to the xref table in an ensembl core database'''


class GeneStableIdAdaptor(sqlgraph.SQLTable):
    '''Provides access to the gene_stable_id table in an ensembl core database'''

class TranscriptStableIdAdaptor(sqlgraph.SQLTable):
    '''Provides access to the transcript_stable_id table in an ensembl core database'''

    
class TranslationStableIdAdaptor(sqlgraph.SQLTable):
    '''Provides access to the translation_stable_id table in an ensembl core database'''

    
class ExonStableIdAdaptor(sqlgraph.SQLTable):
    '''Provides access to the exon_stable_id table in an ensembl core database'''


class SeqregionAdaptor(sqlgraph.SQLTable):
    '''Provides access to the seq_region table in an ensembl core database'''


class DnaAdaptor(sqlgraph.SQLTable):
    '''Provides access to the dna table in an ensembl core database'''

    def __init__(self, name, **kwargs):
        sqlgraph.SQLTable.__init__(self, name, itemSliceClass=seqdb.SeqDBSlice, attrAlias=dict(seq='sequence'), **kwargs)

    
class PeptideArchiveAdaptor(sqlgraph.SQLTable):
    '''Provides access to the peptide_archive table in an ensembl core database'''

   
"""
class GeneArchiveAdaptor(Adaptor):
    '''Provides access to the gene_archive table in an ensembl core database'''

    def __init__(self, dbname, cursor):
        Adaptor.__init__(self, dbname, 'gene_archive', sqlgraph.TupleO, cursor)
"""


class PredictionExonAdaptor(FeatureAdaptor):
    '''Provides access to the prediction_exon table in an ensembl core database'''

   
class PredictionTranscriptAdaptor(FeatureAdaptor):
    '''Provides access to the prediction_transcript table in an ensembl core database'''


    def fetch_by_display_label(self, display_label):
        '''Obtain prediction_transcript objects by a given display_label
 
        >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> ptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
        >>> pts = ptAdaptor.fetch_by_display_label('GENSCAN00000036948')
        >>> for index, pt in enumerate(pts):
        ...     print 'prediction_transcript', index, ':'
        ...     print 'start:', pt.start
        ...     print 'seq_region_start:', pt.seq_region_start
        ...
        prediction_transcript 0 :
        start: 48848
        seq_region_start: 48849
        '''
            
        t = self.select('where display_label = %s', (display_label))
        prediction_transcripts = []
        for row in t:
            prediction_transcript = self[row.id]
            prediction_transcripts.append(prediction_transcript)
        return prediction_transcripts
    

class ExonAdaptor(FeatureAdaptor):
    '''Provides access to the exon table in an ensembl core database

    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> exonAdaptor = coreDBAdaptor.get_adaptor('exon')
    >>> exons = exonAdaptor.fetch_by_stable_id('ENSE00001493538')
    >>> for e in exons:
    ...     print e.id, e.get_stable_id()
    ...
    95160 ENSE00001493538
    >>> slice = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', end=100000, strand = -1)
    >>> exons = exonAdaptor.fetch_all_by_slice(slice)
    >>> for e in exons:
    ...     print e.exon_id, e.seq_region_id, e.sequence.id, e.seq_region_start, e.sequence.start, e.seq_region_strand, e.orientation, e.get_stable_id()
    ...
    13 226034 chr1 4274 -4365 -1 1 ENSE00001367146
    12 226034 chr1 4863 -4901 -1 1 ENSE00001383334
    11 226034 chr1 5659 -5764 -1 1 ENSE00001388009
    10 226034 chr1 5767 -5810 -1 1 ENSE00001375216
    9 226034 chr1 6470 -6608 -1 1 ENSE00001413614
    8 226034 chr1 6611 -6628 -1 1 ENSE00001401257
    7 226034 chr1 6721 -6918 -1 1 ENSE00001400263
    6 226034 chr1 7096 -7227 -1 1 ENSE00001364065
    5 226034 chr1 7465 -7605 -1 1 ENSE00001405391
    4 226034 chr1 7778 -7924 -1 1 ENSE00001386644
    3 226034 chr1 8131 -8229 -1 1 ENSE00001427794
    2 226034 chr1 14600 -14754 -1 1 ENSE00001424900
    1 226034 chr1 19397 -19669 -1 1 ENSE00001388796
    16 226034 chr1 24417 -25037 -1 1 ENSE00001481229
    15 226034 chr1 25140 -25344 -1 1 ENSE00001481230
    14 226034 chr1 25584 -25944 -1 1 ENSE00001481231
    17 226034 chr1 42912 -42930 1 -1 ENSE00001481227
    18 226034 chr1 44693 -44799 1 -1 ENSE00001429474
    281015 226034 chr1 52878 -53750 1 -1 ENSE00001481223
    5176 226034 chr1 58954 -59871 1 -1 ENSE00001248806
    '''

    def fetch_by_stable_id(self, exonStableID):
        
        dbParams = _get_db_parameters(self.name)
        dbAdaptor = _get_DB_adaptor(dbParams[0], dbParams[1], dbParams[2])
        esAdaptor = dbAdaptor.get_adaptor('exon_stable_id')
        t = esAdaptor.select('where stable_id = %s', (exonStableID))
        exons = []
        for row in t:
            exon = self[row.exon_id]
            exons.append(exon)
        return exons
"""   
    def fetch_exons_by_translation(self, transcript_id, start_exon_id, end_exon_id):
        
        cursor = self.cursor
        n = cursor.execute('select rank from %s.exon_transcript where transcript_id = %%s and exon_id = %%s' %(self.db), (transcript_id, start_exon_id))
        t = cursor.fetchall()
        if n != 1:
            raise KeyError('Warning: duplicated!')
        start_rank = t[0][0]
        n = cursor.execute('select rank from %s.exon_transcript where transcript_id = %%s and exon_id = %%s' %(self.db), (transcript_id, end_exon_id))
        t = cursor.fetchall()
        if n != 1:
            raise KeyError('Warning: duplicated!')
        end_rank = t[0][0]
        
        n = cursor.execute('select exon_id from %s.exon_transcript where transcript_id = %%s and rank >= %%s and rank <= %%s' %(self.db), (transcript_id, start_rank, end_rank))
        t = cursor.fetchall()
        exons = []
        if n == 0:
            return exons
        for row in t:
            e = Exon(row[0])
            exons.append(e)
        return exons
"""

class TranscriptAdaptor(FeatureAdaptor):
    '''Provides access to the transcript table in an ensembl core database
    
    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
    >>> transcripts = transcriptAdaptor.fetch_by_stable_id('ENST00000382841')
    >>> for t in transcripts:
    ...     print t.id, t.get_stable_id()
    ...
    15960 ENST00000382841
    '''
   
    def fetch_by_stable_id(self, transcriptStableID):
        
        dbParams = _get_db_parameters(self.name)
        dbAdaptor = _get_DB_adaptor(dbParams[0], dbParams[1], dbParams[2])
        tsAdaptor = dbAdaptor.get_adaptor('transcript_stable_id')
        t = tsAdaptor.select('where stable_id = %s', (transcriptStableID))
        transcripts = []
        for row in t:
            transcript = self[row.transcript_id]
            transcripts.append(transcript)
        return transcripts



class GeneAdaptor(FeatureAdaptor):
    '''Provides access to the gene table in an ensembl core database

    >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    >>> geneAdaptor = coreDBAdaptor.get_adaptor('gene')
    >>> genes = geneAdaptor.fetch_by_stable_id('ENSG00000120645')
    >>> for g in genes:
    ...     print g.id, g.get_stable_id()
    ...
    8946 ENSG00000120645
    '''


    def fetch_by_stable_id(self, geneStableID):
        
        dbParams = _get_db_parameters(self.name)
        dbAdaptor = _get_DB_adaptor(dbParams[0], dbParams[1], dbParams[2])
        gsAdaptor = dbAdaptor.get_adaptor('gene_stable_id')
        t = gsAdaptor.select('where stable_id = %s', (geneStableID))
        genes = []
        for row in t:
            gene = self[row.gene_id]
            genes.append(gene)
        return genes

    """
    def fetch_genes_by_externalRef(self, external_ref_label):
        '''Return all the genes that are associated with the given external reference or an empty set.'''
        
        '''
        cursor = self.cursor
        #print('select xref_id from %s.xref where display_label = %s' %(self.db, external_ref_label)) 
        
        n = cursor.execute('select xref_id from %s.xref where display_label = %%s' %(self.db), (external_ref_label))
        t = cursor.fetchall()
        '''
        # Retrieve xref_id(s) that are associated with the given external reference
        
        driver = getDriver('ensembldb.ensembl.org', 'anonymous', 'homo_sapiens_core_47_36i')
        xref_adaptor = driver.getAdaptor('xref')
        t = xref_adaptor.tbobj.select('where display_label = %s', (external_ref_label), None, 't1.xref_id')
        
        # Retrieve gene_id(s) that are associated with the returned xref_id(s)
        geneIDs = []
        for row in t:
            xref_id = row.xref_id
            #print xref_id
            gene_rows = self.tbobj.select('where display_xref_id = %s', (xref_id), None, 't1.gene_id')
            for gid_row in gene_rows:
                geneIDs.append(gid_row.gene_id)
        # Create genes based on the retrieved gene_ids
        genes = []
        for gid in geneIDs:
            g = Gene(gid)
            genes.append(g)
        return genes
   """


class CoordSystemAdaptor(sqlgraph.SQLTable):
    'Provide access to the coord_system table in an ensembl core database'


class CoreDBAdaptor(object):
    'Provide access to a core database in the ensembl database system'

    # static attributes
    dbType = 'core'
    TBAdaptorClass = {'meta_coord': MetaCoordAdaptor, 'dna': DnaAdaptor, 'exon': ExonAdaptor, 'gene': GeneAdaptor, 'transcript': TranscriptAdaptor, 'seq_region': SeqregionAdaptor, 'translation': TranslationAdaptor, 'gene_stable_id': GeneStableIdAdaptor, 'transcript_stable_id': TranscriptStableIdAdaptor, 'translation_stable_id': TranslationStableIdAdaptor, 'exon_stable_id': ExonStableIdAdaptor, 'xref': XrefAdaptor, 'prediction_transcript': PredictionTranscriptAdaptor, 'meta_coord': MetaCoordAdaptor, 'prediction_exon': PredictionExonAdaptor, 'peptide_archive': PeptideArchiveAdaptor, 'coord_system': CoordSystemAdaptor}
    RowClass = {'meta_coord': MetaCoord, 'dna': Dna, 'exon': Exon, 'gene': Gene, 'transcript': Transcript, 'seq_region': Seqregion, 'translation': Translation, 'gene_stable_id': GeneStableID, 'transcript_stable_id': TranscriptStableID, 'translation_stable_id': TranslationStableID, 'exon_stable_id': ExonStableID, 'xref': Xref, 'prediction_transcript': PredictionTranscript, 'prediction_exon': PredictionExon, 'peptide_archive': PeptideArchive, 'coord_system': CoordSystem}
    AnnotClass = {'exon': ExonAnnot, 'transcript': TranscriptAnnot, 'gene': GeneAnnot}
    #MapperClass = {'prediction_transcript_prediction_exon': PtranscriptToPexon, 'gene_transcript': GeneToTranscript, 'transcript_translation': TranscriptToTranslation}
    #Genomes = {'homo_sapiens_47_36i': 'HUMAN.hg18'}

    #import pygr.Data
    RelationClass = {'ManyToMany': pygr.Data.ManyToManyRelation, 'OneToMany': pygr.Data.OneToManyRelation, 'OneToOne': pygr.Data.OneToOneRelation}
    #GraphClass = {'homo_sapiens_core_47_36i.exon_transcript': pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.exon_transcript}
    #SchemaClass = {'homo_sapiens_core_47_36i.exon_transcript': pygr.Data.schema.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.exon_transcript}

    #GraphClass = {'homo_sapiens_core_47_36i': pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i}
    #SchemaClass = {'homo_sapiens_core_47_36i': pygr.Data.schema.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i}

    def __init__(self, conn, dbSpecies, dbVersion):
        
        # instance attributes
        self.conn = conn
        self.dbSpecies = dbSpecies
        self.dbVersion = dbVersion
        self.dbName = self.dbSpecies+'_'+CoreDBAdaptor.dbType+'_'+self.dbVersion
        self.genomeName = 'HUMAN.hg18'
        self.graphClass = pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i
        self.schemaClass = pygr.Data.schema.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i
    
    def get_feature(self, featureTB):
        'Obtain an annotation DB for a core database table'
        # get a SQLTable object
        name = self.dbName + '.' + featureTB
        #tbAdaptor = sqlgraph.SQLTable(name, itemClass=AnnotBaseModel, serverInfo=self.conn)
        tbAdaptor = sqlgraph.SQLTable(name, itemClass=seqregion.EnsemblRow, serverInfo=self.conn)
        annoDB = self._get_annotationDB(featureTB, tbAdaptor)
        return annoDB
            
    def get_adaptor(self, tbname): 
        'Obtain a particular table adaptor of a core database table'

        name = self.dbName + '.' + tbname
        # Get the tbobj from pygr.Data
        resource_id = 'Bio.Annotation.Ensembl.' + name + '.' + 'sqltable'
        #print self.resource_id
        tbAdaptor = _get_resource(resource_id)
        if tbAdaptor == None:
            # Create a pickleable SQLTable adaptor
            adaptorClass = CoreDBAdaptor.TBAdaptorClass[tbname]
            rowClass = CoreDBAdaptor.RowClass[tbname]
            tbAdaptor = adaptorClass(name, itemClass=rowClass, serverInfo=self.conn)
           
            # Save self.tbobj to pygr.Data
            _save_resource(resource_id, tbAdaptor, 'sqltable')
        return tbAdaptor


    def _create_featureMapper(self, sourceTBName, targetTBName):

        tableResourceID = 'Bio.Annotation.Ensembl'
        sourceResID = tableResourceID + '.' + self.dbName + '.' + sourceTBName + '.' + 'sqltable'
        #print 'sourceResID: ', sourceResID
        targetResID = tableResourceID + '.' + self.dbName + '.' + targetTBName + '.' + 'sqltable'
        #print 'targetResID: ', targetResID
        sourceTB = _get_resource(sourceResID)
        if sourceTB == None:
            sourceTB = self.get_adaptor(sourceTBName)
            _save_resource(sourceResID, sourceTB, 'sqltable')
        targetTB = _get_resource(targetResID)
        if targetTB == None:
            targetTB = self.get_adaptor(targetTBName)
            _save_resource(targetResID, targetTB, 'sqltable')
        mapperName = sourceTBName + '_' + targetTBName
        mapperClass = CoreDBAdaptor.MapperClass[mapperName] 
        myMapper = mapperClass(sourceTB, targetTB)
        mapperData = [myMapper, sourceTB, targetTB]
        return mapperData

    def _fetch_featureMapper(self, sourceTBName, targetTBName):
        mapperName = sourceTBName + '_' + targetTBName
        mapperID = 'Bio.Annotation.Ensembl.' + self.dbName + '.' + mapperName
        mapper = _get_resource(mapperID)
        if mapper == None:
            mapperData = self._create_featureMapper(sourceTBName, targetTBName)
            mapper = mapperData[0]
            sourceTB = mapperData[1]
            targetTB = mapperData[2]
            #_save_mapper(mapperID, mapper, sourceTB, targetTB)
        return mapper

    def _create_graph(self, source, target, relationTBName):

        sourceTB = source.sliceDB
        sourceTBName = sourceTB.name
        sourceTBNameList = sourceTBName.split('.')
        sourceName = sourceTBNameList[1]
        targetTB = target.sliceDB
        targetTBName = targetTB.name
        targetTBNameList = targetTBName.split('.')
        targetName = targetTBNameList[1]
        spliceName = self.dbName + '.' + relationTBName
        # test
        #print spliceName
        
        sourceID = sourceName + '_id'
        targetID = targetName + '_id'
        # test
        #print sourceID
        #print targetID
        # create a forward mapping from sourceDB -> targetDB
        myGraph = sqlgraph.SQLGraph(spliceName,serverInfo=self.conn, sourceDB=source, targetDB=target, attrAlias=dict(source_id=sourceID, target_id=targetID))
        return myGraph
        
    def _create_featureGraph(self, sourceTBName, targetTBName):
        tableResourceID = 'Bio.Annotation.Ensembl'
    
        sourceResID = tableResourceID + '.' + self.dbName + '.' + sourceTBName + '.' + 'sqltable'
        #print 'sourceResID: ', sourceResID
        targetResID = tableResourceID + '.' + self.dbName + '.' + targetTBName + '.' + 'sqltable'
        #print 'targetResID: ', targetResID
        sourceTB = _get_resource(sourceResID)
        if sourceTB == None:
            sourceTB = self.get_adaptor(sourceTBName)
            _save_resource(sourceResID, sourceTB, 'sqltable')
        targetTB = _get_resource(targetResID)
        if targetTB == None:
            targetTB = self.get_adaptor(targetTBName)
            _save_resource(targetResID, targetTB, 'sqltable')
        graphName = sourceTBName + '_' + targetTBName
        name = self.dbName + '.' + graphName 
        sourceAnnoDB = self._get_annotationDB(sourceTBName, sourceTB)
        targetAnnoDB = self._get_annotationDB(targetTBName, targetTB)
        # create a forward mapping from sourceDB -> targetDB
        myGraph = sqlgraph.SQLGraph(name,serverInfo=self.conn, sourceDB=sourceAnnoDB, targetDB=targetAnnoDB, attrAlias=dict(source_id='exon_id', target_id='transcript_id')) 
        graphData = [myGraph, sourceAnnoDB, targetAnnoDB]
        return graphData
    
    def _save_graph(self, graph, source, target, relation, sourceAttr, targetAttr):
        '''Save the graph to pygr.Data'''

        sourceTB = source.sliceDB
        sourceTBName = sourceTB.name
        sourceTBNameList = sourceTBName.split('.')
        sourceName = sourceTBNameList[1]
        targetTB = target.sliceDB
        targetTBName = targetTB.name
        targetTBNameList = targetTBName.split('.')
        targetName = targetTBNameList[1]
        graphKey = self.dbName + '.' + sourceName + '_' + targetName
        relationKey = relation 
        doc = 'ensembl annotation graph ' + sourceName + ' -> ' + targetName + ' (' + self.dbName + ')'
        graph.__doc__ = doc
        #import pygr.Data
        #graphClass = self.GraphClass[graphKey]
        graphClass = self.graphClass
        #schemaClass = self.SchemaClass[graphKey]
        schemaClass = self.schemaClass
        relationClass = self.RelationClass[relationKey]
        #relationClass = self.RelationClass[relationKey]
        schemaValue = relationClass(source, target, bindAttrs=(sourceAttr, targetAttr))
        attribute = sourceName + '_' + targetName
        #print attribute
        graphClass.__setattr__(attribute, graph)
        schemaClass.__setattr__(attribute, schemaValue)
        pygr.Data.save()
    
    def _save_featureGraph(self, graph, sourceTBName, sourceAnnoDB, targetTBName, targetAnnoDB, relation, sourceAttr, targetAttr):
        '''Save the graph to pygr.Data'''

        graphKey = self.dbName + '.' + sourceTBName + '_' + targetTBName
        relationKey = relation 
        doc = 'ensembl annotation graph ' + sourceTBName + ' -> ' + targetTBName + ' (' + self.dbName + ')'
        graph.__doc__ = doc
        import pygr.Data
        graphClass = self.GraphClass[graphKey]
        schemaClass = self.SchemaClass[graphKey]
        relationClass = self.RelationClass[relationKey]
        schemaValue = relationClass(sourceAnnoDB, targetAnnoDB, bindAttrs=(sourceAttr, targetAttr))
        attribute = sourceTBName + '_' + targetTBName
        graphClass.__setattr__(attribute, graph)
        schemaClass.__setattr__(attribute, schemaValue)
        pygr.Data.save()

        '''
    def _fetch_featureGraph(self, sourceTBName, targetTBName):
        'Obtain a feature -> feature graph object and its associated sourceAnnotationDB and targetAnnotationDB'

        graphName = sourceTBName + '_' + targetTBName
        #graphName = targetTBName + '_' + targetTBName
        graphID = 'Bio.Annotation.Ensembl.' + self.dbName + '.' + graphName
        #print graphID
        graph = _get_resource(graphID)
        if graph == None:
            # print 'None'
            # create a SQLGraph object
            graphData = self._create_featureGraph(sourceTBName, targetTBName)
            graph = graphData[0]
            sourceDB = graphData[1]
            targetDB = graphData[2]
            
            #_save_graph(graphID, graph, sourceDB, targetDB)
        
        else:
            print 'graph already exists'
            tableResourceID = 'Bio.Annotation.Ensembl'
    
            sourceResID = tableResourceID + '.' + self.dbName + '.' + sourceTBName + '.' + 'sqltable'
            #print 'sourceResID: ', sourceResID
            targetResID = tableResourceID + '.' + self.dbName + '.' + targetTBName + '.' + 'sqltable'
            #print 'targetResID: ', targetResID
            sourceTB = _get_resource(sourceResID)
            if sourceTB == None:
                sourceTB = self.get_adaptor(sourceTBName)
                _save_resource(sourceResID, sourceTB, 'sqltable')
            targetTB = _get_resource(targetResID)
            if targetTB == None:
                targetTB = self.get_adaptor(targetTBName)
                _save_resource(targetResID, targetTB, 'sqltable')
            sourceAnnoDB = self._get_annotationDB(sourceTBName, sourceTB)
            targetAnnoDB = self._get_annotationDB(targetTBName, targetTB)
            graphData = [graph, sourceAnnoDB, targetAnnoDB]
         
        return graphData
        '''            

    def _create_seqregion(self):
        
        #srDBAdaptor = self.get_adaptor('seq_region')
        seq_region_tbname = self.dbName + '.seq_region'
        srDBAdaptor = sqlgraph.SQLTable(seq_region_tbname, itemClass=sqlgraph.TupleO, serverInfo=self.conn)
        dna_tbname = self.dbName + '.dna'
        #dnaAdaptor = self.get_adaptor('dna')
        dnaAdaptor = sqlgraph.SQLTable(dna_tbname, itemSliceClass=seqdb.SeqDBSlice, attrAlias=dict(seq='sequence'), itemClass=EnsemblDNA, serverInfo=self.conn)  
        #genomeKey = self.dbSpecies + '_' + self.dbVersion
        #print genomeKey
        #genomeName = CoreDBAdaptor.Genomes[genomeKey]
        genomeResourceID = 'Bio.Seq.Genome.' + self.genomeName
        #print genomeResourceID
        genome = _get_resource(genomeResourceID)
        #chr1seq = genome['chr1']
        #print repr(chr1seq)
        # create a SeqRegion object
        sr = SeqRegion(srDBAdaptor, {17:genome, 4:dnaAdaptor}, {17:'chr', 4:None})
        return sr

    def _get_seqregion(self):
        # get the seqregion object from pygr.Data
        srID = 'Bio.Seq.Ensembl.' + self.dbName + '.seq'
        seqRegion = _get_resource(srID)
        if seqRegion == None:
            # create a seqRegionDB
            seqRegion = self._create_seqregion()
            # save the seqRegionDB
            _save_resource(srID, seqRegion)
        return seqRegion

    def _create_assemblyMapper(self):
        sr = self._get_seqregion()
        amap = AssemblyMapper(sr, 4, 17)
        return amap

    def _get_assemblyMapper(self):
        # get the assemblyMapper from pygr.Data
        amapID = 'Bio.MSA.Ensembl.' + self.dbName + '.assembly'
        amap = _get_resource(amapID)
        if amap == None:
            # create an assemblyMapper
            amap = self._create_assemblyMapper()
            # save the assemblyMapper
            _save_resource(amapID, amap)
        return amap
            
    def _create_annoDB(self, tbAdaptor):
        sr = self._get_seqregion()
        #from pygr.seqdb import AnnotationDB
        tbName = tbAdaptor.name.split('.')[1]
        annoDB = seqdb.AnnotationDB(tbAdaptor, sr, itemClass=CoreDBAdaptor.AnnotClass[tbName], sliceAttrDict= dict(id='seq_region_id', stop='seq_region_end', orientation='seq_region_strand'))
        return annoDB
        
    def _get_annotationDB(self, tbname, tbAdaptor):
        # get the annotation DB from pygr.Data
        annoDBID = 'Bio.Annotation.Ensembl.' + self.dbName + '.' + tbname
        #print annoDBID
        annoDB = _get_resource(annoDBID)
        if annoDB == None:
            # create a annoDB
            annoDB =self._create_annoDB(tbAdaptor)
            # save the annoDB
            _save_resource(annoDBID, annoDB)
        return annoDB

    def _create_seqFeatureMapper(self, tbname, tbAdaptor):
        annoDB = self._get_annotationDB(tbname, tbAdaptor)
        sr = self._get_seqregion()
        mapper = EnsemblMapper(annoDB, sr)
        mapperData = [mapper, annoDB, sr]
        return mapperData
        
        
    def _get_seqFeatureMapper(self, tbName, tbAdaptor):
        # get the seqFeatureMapper from pygr.Data
        seqFeatureMapperID = 'Bio.Annotation.Ensembl.' + self.dbName + '.' + 'slice' + 'To' + tbName + 's'
        #print seqFeatureMapperID
        seqFeatureMapper = _get_resource(seqFeatureMapperID)
        if seqFeatureMapper == None:
            # create a seqFeatureMapper
            mapperData = self._create_seqFeatureMapper(tbName, tbAdaptor)
            seqFeatureMapper = mapperData[0]
            annoDB = mapperData[1]
            sr = mapperData[2]
        
            # save the seqFeatureMappper
            #_save_mapper(mapperID, seqFeatureMapper, sr, annoDB)
            
        return seqFeatureMapper

    
    def fetch_slice_by_region(self, coordSystemName, seqregionName, start=None, end=None, strand=None):
        '''Obtain the DNA sequence of a particular genomic region (defined by
chromosome, chrName, start, end, strand) or a particular contig region (defined by contig, contigName, start, end, strand).
        Note: the start and end are based on ensembl 1-offset coordinate system.

        >>> serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
        >>> coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
        >>> s = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', start=1000001, end=1000010, strand=-1)
        >>> print str(s)
        GCAGCCACGT
        >>> s = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', 1000001, 1000010)
        >>> print str(s)
        ACGTGGCTGC
        >>> s = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', end=10)
        >>> print str(s)
        taaccctaac
        >>> s = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', start=1, end=10)
        >>> print str(s)
        taaccctaac
        >>> s = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', start=247249710)
        >>> print str(s)
        NNNNNNNNNN
        >>> s = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', start=247249710, end=247249720)
        >>> print str(s)
        NNNNNNNNNN
        >>> s = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', strand=-1)
        >>> print repr(s)
        -chr1[0:247249719]
        >>> s = coreDBAdaptor.fetch_slice_by_region('contig', 'AADC01095577.1.1.41877') 
        >>> print repr(s)
        143909[0:41877]
        >>> print len(s)
        41877
        >>> s = coreDBAdaptor.fetch_slice_by_region('contig', 'AADC01095577.1.1.41877', end=10)
        >>> print str(s)
        CACCCTGCCC
        '''
        
        # get the seqregion DB
        sr = self._get_seqregion()
        # retrieve the required seq_region_id from the coord_system table and the seq_region table by the coordSystemName and seqregionName
        coordSystemAdaptor = self.get_adaptor('coord_system')
        t1 = coordSystemAdaptor.select('where name = %s', (coordSystemName))
       
        coordSystemID = t1.next().coord_system_id
        if coordSystemID == 101:
            coordSystemID = 17 # coord_system_id for 'chromosome'
        #print coordSystemID
        t2 = sr.seqRegionDB.select('where name = %s and coord_system_id = %s' , (seqregionName, coordSystemID))
        seq_regionID = t2.next().seq_region_id
        # get the entire sequence defined by the seq_regionID
        slice = sr[seq_regionID]
        # convert an ensembl start coordinate to a Python zero-off coordinate
        if start is not None:
            python_start = start - 1
            if end is None:
                # obtain the required interval
                slice = slice[python_start:]
            else:
                slice = slice[python_start:end]
        else:
            if end is not None:
                slice = slice[:end]
        # If the required interval is on the reverse strand, then return the reverse and complimented sequence interval.
        if strand == -1:
            slice = -slice
        return slice
    
def _test():
    import doctest
    doctest.testmod()
    
if __name__ == '__main__': # example code
    
    #_test()
    
    serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    exonDB = coreDBAdaptor.get_feature('exon')
    exon = exonDB[1]
    print str(exon.sequence)
    #slice = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', end=100000, strand = -1)
    #print slice.id
    #slice = coreDBAdaptor.fetch_slice_by_region('chromosome', '18', start=31129340, end=31190706)
    '''
    slice = coreDBAdaptor.fetch_slice_by_region('contig', 'AC116447.10.1.105108') 
    ptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
    pts = ptAdaptor.fetch_all_by_slice(slice)
    for pt in pts:
        print repr(pt), pt.id, pt.display_label, pt.sequence.id, pt.sequence.start, pt.seq_region_start, pt.sequence.stop, pt.seq_region_end, pt.seq_region_strand, pt.orientation
    
    #s = coreDBAdaptor.fetch_slice_by_region('chromosome', 1, 1000001, 1000010)
    #s = coreDBAdaptor.fetch_slice_by_region('chromosome', 1)
    #s = coreDBAdaptor.fetch_slice_by_region('chromosome', 1, start=1000001, end=1000010, strand=-1)
    s = coreDBAdaptor.fetch_slice_by_region('contig', 'AADC01095577.1.1.41877') 
    print repr(s) 
    print len(s)
    s = coreDBAdaptor.fetch_slice_by_region('contig', 'AADC01095577.1.1.41877', end=10)
    print str(s)
    #seqregionAdaptor = coreDBAdaptor.get_adaptor('seq_region')
    #dna = sqlgraph.SQLTable('homo_sapiens_core_47_36i.dna', serverInfo=conn, 
    #                        itemClass=EnsemblDNA,
    #                        itemSliceClass=seqdb.SeqDBSlice,
    #                        attrAlias=dict(seq='sequence'))
    #hg18 = _get_resource('Bio.Seq.Genome.HUMAN.hg18')
    #srdb = SeqRegion(seq_region, {17:hg18, 4:dna}, {17:'chr',4:None})
    
    
    peAdaptor = coreDBAdaptor.get_adaptor('prediction_exon')
    predictionExons = peAdaptor.fetch_all_by_slice(slice)
    for pe in predictionExons:
        print pe.id
    
    geneAdaptor = coreDBAdaptor.get_adaptor('gene')
    genes = geneAdaptor.fetch_all_by_slice(slice)
    for g in genes:
        print g.gene_id, g.seq_region_id, g.sequence.id, g.seq_region_strand, g.orientation, g.get_stable_id()
    
    transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
    transcripts = transcriptAdaptor.fetch_all_by_slice(slice)
    for t in transcripts:
        print t.transcript_id, t.seq_region_id, t.sequence.id, t.seq_region_strand, t.orientation, t.get_gene().gene_id, t.get_translation().translation_id
    
    exonAdaptor = coreDBAdaptor.get_adaptor('exon')
    #exonAnnoDB = coreDBAdaptor._get_annotationDB('exon', exonAdaptor)
    exon = exonAdaptor[73777]
    exonSeq = exon.get_sequence()
    print str(exonSeq), exonSeq.orientation, repr(exonSeq), len(exonSeq)
    exonSlice = exon.get_sequence(10)
    print str(exonSlice), exonSlice.orientation, repr(exonSlice), len(exonSlice)
    #annoExon = exonAnnoDB[73777]
    #print repr(annoExon), repr(annoExon.sequence), annoExon.sequence.orientation, str(annoExon.sequence)
    #slice = coreDBAdaptor.fetch_slice_by_region('chromosome', '10', start=444866, end=444957)
    #print str(slice), slice.orientation, len(slice)
    
    ptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
    pt = ptAdaptor[36948]
    ptSeq = pt.get_sequence()
    print repr(ptSeq)
    
    exons = exonAdaptor.fetch_all_by_slice(slice)
    #print exons
    for e in exons:
        print e.exon_id, e.seq_region_id, e.sequence.id, e.seq_region_start, e.sequence.start, e.seq_region_strand, e.orientation, e.get_stable_id()
    
    #transcriptToExons = TranscriptToExon(transcriptAdaptor, exonAdaptor)
    transcript = transcriptAdaptor[1]
    #exons = transcriptToExons[transcript]
    exons = transcript.get_all_exons()
    print 'Transcript', transcript.id, 'has', len(exons), 'exons:'
    print 'id', 'start'
    for e in exons:
        print e.id, e.start
    exon = exonAdaptor[1]
    #transcripts = (~transcriptToExons)[exon]
    transcripts = exon.get_all_transcripts()
    print 'Exon', exon.id, 'belongs to', len(transcripts), 'transcript:'
    print 'id', 'start'
    for t in transcripts:
        print t.id, t.start
    
    translationAdaptor = coreDBAdaptor.get_adaptor('translation')
    translation = translationAdaptor[10]
    print translation.transcript_id
    print translation.start_exon_id
    print translation.end_exon_id
    translations = translationAdaptor.fetch_by_stable_id('ENSP00000372292')
    for tl in translations:
        print tl.id, tl.get_stable_id()

    geneAdaptor = coreDBAdaptor.get_adaptor('gene')
    genes = geneAdaptor.fetch_by_stable_id('ENSG00000120645')
    for g in genes:
        print g.id, g.get_stable_id()

    transcriptAdaptor = coreDBAdaptor.get_adaptor('transcript')
    transcripts = transcriptAdaptor.fetch_by_stable_id('ENST00000382841')
    for t in transcripts:
        print t.id, t.get_stable_id()

    exonAdaptor = coreDBAdaptor.get_adaptor('exon')
    exons = exonAdaptor.fetch_by_stable_id('ENSE00001493538')
    for e in exons:
        print e.id, e.get_stable_id()
    #translation.getAttributes()
    #exons = translation.getExons()
    #print 'id', 'start'
    #for e in exons:
    #    print e.id, e.start
    
    
    ptAdaptor = coreDBAdaptor.get_adaptor('prediction_transcript')
    pts = ptAdaptor.fetch_by_display_label('GENSCAN00000036948')
    
    for index, pt in enumerate(pts):
        print 'prediction_transcript', index, ':'
        pt.getAttributes()
        print 'start:', pt.start
        print 'seq_region_id:', pt.seq_region_id
        print 'seq_region_start:', pt.seq_region_start
        print 'seq_region_end:', pt.seq_region_end
        print 'seq_region_strand:', pt.seq_region_strand
        print 'analysis_id:', pt.analysis_id
        print 'display_label:', pt.display_label
    '''
    '''
    exon_adaptor = driver.getAdaptor('exon')
    #exons = exon_adaptor.fetch_exons_by_seqregion(1, 6023217, 6023986, 1, driver)
    #exons = exon_adaptor.fetch_exons_by_seqregion(10, 444866, 444957, -1, driver)
    exons = exon_adaptor.fetch_exons_by_seqregion(1, 1, 100000, 1, driver)
    for index, e in enumerate(exons):
       print '\nexon', index 
       e.getAttributes()
    
    print '\ngene_adaptor.fetch_genes_by_externalRef():'
    gene_adaptor = driver.getAdaptor('gene')
    genes = gene_adaptor.fetch_genes_by_externalRef('IQSEC3')
    #genes = gene_adaptor.fetch_genes_by_externalRef('GO:0000228')
    if len(genes) == 0:
        print '\nNo gene is associated with the given external reference label.'
    else:
        for index, g in enumerate(genes):
            print '\ngene', index, ':'
            g.getAttributes()
       
    genes = gene_adaptor.fetch_genes_by_seqregion(1, 4274, 19669, -1, driver)
    for index, g in enumerate(genes):
        print '\ngene', index
        g.getAttributes()
        g.getSequence('gene')
    
    transcript_adaptor = driver.getAdaptor('transcript')
    
    transcripts = transcript_adaptor.fetch_transcripts_by_seqregion(1, 4274, 19669, -1, driver) 
    for index, t in enumerate(transcripts):
        print '\ntranscript', index
        t.getAttributes()
        #t.getSequence('transcript')
    
    print '\ntranscript_adaptor.fetch_transcripts_by_geneID(gene_id):'
    transcripts = transcript_adaptor.fetch_transcripts_by_geneID(8946)
    if len(transcripts) == 0:
        print '\nNo transcript identified for this gene.'
    else:
        for index, t in enumerate(transcripts):
            print '\ntranscript ', index, ':'
            t.getAttributes()
            print 'length: ', len(t.getSequence('transcript'))
    
    
    #s = driver.fetch_sequence_by_region(1, 4274, 19669, 1)
    s = driver.fetch_sequence_by_region(10, 444866, 444957, 1)
    print "\nLength of the sequence: ", len(s)
    print "The sequence: ", str(s)
   
   
    gene_stableID_adaptor = driver.getAdaptor('gene_stable_id')
    print "\ntest gene_stableID_adaptor.fetch_by_stable_id('ENSG00000215911')"
    #genes = gene_stableID_adaptor.fetch_by_stable_id('ENSG00000215911')
    #_retrieve_units_tester(genes, 'gene')
    _fetch_by_stableID_tester(gene_stableID_adaptor, 'ENSG00000215911')

    transcript_stableID_adaptor = driver.getAdaptor('transcript_stable_id')
    print "\ntest transcript_stableID_adaptor.fetch_by_stable_id('ENST00000382841')"
    _fetch_by_stableID_tester(transcript_stableID_adaptor, 'ENST00000382841')

    translation_stableID_adaptor = driver.getAdaptor('translation_stable_id')
    print "\ntest translation_stableID_adaptor.fetch_by_stable_id('ENSP00000317958')"
    _fetch_by_stableID_tester(translation_stableID_adaptor, 'ENSP00000317958')

    exon_stableID_adaptor = driver.getAdaptor('exon_stable_id')
    print "\ntest exon_stableID_adaptor.fetch_by_stable_id('ENSE00001493538')"
    _fetch_by_stableID_tester(exon_stableID_adaptor, 'ENSE00001493538')
    
    
    prediction_transcript_adaptor = driver.getAdaptor('prediction_transcript')
    prediction_transcripts = prediction_transcript_adaptor.fetch_by_display_label('GENSCAN00000036948')
    for index, pt in enumerate(prediction_transcripts):
        print 'prediction_transcript', index, ':'
        pt.getAttributes()
        print 'start:', pt.rowobj.start
        print 'seq_region_id:', pt.getSeqregionID()
        print 'seq_region_start:', pt.getSeqregionStart()
        print 'seq_region_end:', pt.getSeqregionEnd()
        print 'seq_region_strand:', pt.getOrientation()
        print 'analysis_id:', pt.getAnalysisID()
        print 'display_label:', pt.getDisplayLabel()
    '''   

