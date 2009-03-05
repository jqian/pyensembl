
import sys, os, string

seqDir = '/result/pygr_megatest/seq_data' # SEQDB.BLASTDB
msaDir = '/result/pygr_megatest/maf_test' # PRE BUILT NLMSA

## msaDir CONTAINS PRE-BUILT NLMSA
## seqDir CONTAINS GENOME ASSEMBLIES AND THEIR SEQDB FILES
## TEST INPUT/OUPTUT FOR COMPARISON, THESE FILES SHOULD BE IN THIS DIRECTORY
##        exonAnnotFileName = 'Annotation_ConservedElement_Exons_dm2.txt'
##        intronAnnotFileName = 'Annotation_ConservedElement_Introns_dm2.txt'
##        stopAnnotFileName = 'Annotation_ConservedElement_Stop_dm2.txt'
## testDir = os.path.join('/usr/tmp/deepreds', 'TEST_' + ''.join(tmpList)) SHOULD BE DELETED IF YOU WANT TO RUN IN '.'

# DIRECTIONARY FOR DOC STRING OF SEQDB
docStringDict = {
    'anoGam1':'A. gambiae Genome (February 2003)',
    'apiMel2':'A. mellifera Genome (January 2005)',
    'dm2':'D. melanogaster Genome (April 2004)',
    'dp4':'D. pseudoobscura Genome (February 2006)',
    'droAna3':'D. ananassae Genome (February 2006)',
    'droEre2':'D. erecta Genome (February 2006)',
    'droGri2':'D. grimshawi Genome (February 2006)',
    'droMoj3':'D. mojavensis Genome (February 2006)',
    'droPer1':'D. persimilis Genome (October 2005)',
    'droSec1':'D. sechellia Genome (October 2005)',
    'droSim1':'D. simulans Genome (April 2005)',
    'droVir3':'D. virilis Genome (February 2006)',
    'droWil1':'D. willistoni Genome (February 2006)',
    'droYak2':'D. yakuba Genome (November 2005)',
    'triCas2':'T. castaneum Genome (September 2005)',
    }

# GENOME ASSEMBLY LIST FOR DM2 MULTIZ15WAY
msaSpeciesList = ['anoGam1', 'apiMel2', 'dm2', 'dp4', 'droAna3', 'droEre2', 'droGri2', 'droMoj3', \
                  'droPer1', 'droSec1', 'droSim1', 'droVir3', 'droWil1', 'droYak2', 'triCas2']

class PygrBuildNLMSAMegabase(object):
    'restrict megatest to an initially empty directory, need large space to perform'
    def __init__(self, testDir = None):
        import random
        tmpList = [c for c in 'PygrBuildNLMSAMegabase']
        random.shuffle(tmpList)
        testDir = os.path.join('/usr/tmp/deepreds', 'TEST_' + ''.join(tmpList)) # FOR TEST, SHOULD BE DELETED
        if testDir is None: testDir = 'TEST_' + ''.join(tmpList) # NOT SPECIFIED, USE CURRENT DIRECTORY
        try:
            os.mkdir(testDir)
            testDir = os.path.realpath(testDir)
        except:
            raise IOError
        self.path = testDir
        try:
            tmpFileName = os.path.join(testDir, 'DELETE_THIS_TEMP_FILE')
            open(tmpFileName, 'w').write('A'*1024*1024) # WRITE 1MB FILE FOR TESTING
        except:
            raise IOError
        os.environ['PYGRDATAPATH'] = self.path
        import pygr.Data
        from pygr import seqdb
        for orgstr in msaSpeciesList:
            genome = seqdb.BlastDB(os.path.join(seqDir, orgstr))
            genome.__doc__ = docStringDict[orgstr]
            pygr.Data.getResource.addResource('TEST.Seq.Genome.' + orgstr, genome)
        pygr.Data.save()
    def copyFile(self, filename): # COPY A FILE INTO TEST DIRECTORY
        newname = os.path.join(self.path, os.path.basename(filename))
        open(newname, 'w').write(open(filename, 'r').read())
        return newname
    def teardown(self):
        'delete the temporary directory and files'
        for dirpath, subdirs, files in os.walk(self.path, topdown = False): # SHOULD BE DELETED BOTTOM-UP FASHION
            # THIS PART MAY NOT WORK IN NFS MOUNTED DIRECTORY DUE TO .nfsXXXXXXXXX CREATION
            # IN NFS MOUNTED DIRECTORY, IT CANNOT BE DELETED UNTIL CLOSING PYGRDATA
            for filename in files:
                os.remove(os.path.join(dirpath, filename))
            os.rmdir(dirpath)

class Build_Test(PygrBuildNLMSAMegabase):
    def seqdb_test(self): # CHECK PYGR.DATA CONTENTS
        os.environ['PYGRDATAPATH'] = self.path
        import pygr.Data
        l = pygr.Data.dir('TEST')
        preList = ['TEST.Seq.Genome.' + orgstr for orgstr in msaSpeciesList]
        assert l == preList
    def collectionannot_test(self): # BUILD ANNOTATION DB FROM FILE
        os.environ['PYGRDATAPATH'] = self.path
        import pygr.Data
        from pygr import seqdb, cnestedlist, sqlgraph
        dm2 = pygr.Data.getResource('TEST.Seq.Genome.dm2')
        # BUILD ANNOTATION DATABASE FOR REFSEQ EXONS
        exon_slices = pygr.Data.Collection(filename = os.path.join(self.path, 'refGene_exonAnnot_dm2.cdb'), \
            intKeys = True, mode = 'c', writeback = False) # ONLY C
        exon_db = seqdb.AnnotationDB(exon_slices, dm2,
                               sliceAttrDict = dict(id = 0, exon_id = 1, orientation = 2,
                                                  gene_id = 3, start = 4, stop = 5))
        msa = cnestedlist.NLMSA(os.path.join(self.path, 'refGene_exonAnnot_dm2'), 'w', \
            pairwiseMode = True, bidirectional = False)
        #for lines in open('refGene_exonAnnot_dm2.txt', 'r').xreadlines():
        for lines in open('refGene_exonAnnot_chrYh_dm2.txt', 'r').xreadlines():
            row = [x for x in lines.split('\t')] # CONVERT TO LIST SO MUTABLE
            row[1] = int(row[1]) # CONVERT FROM STRING TO INTEGER
            exon_slices[row[1]] = row
            exon = exon_db[row[1]] # GET THE ANNOTATION OBJECT FOR THIS EXON
            msa.addAnnotation(exon) # SAVE IT TO GENOME MAPPING
        exon_db.clear_cache() # not really necessary; cache should autoGC
        exon_slices.close() # SHELVE SHOULD BE EXPLICITLY CLOSED IN ORDER TO SAVE CURRENT CONTENTS
        msa.build() # FINALIZE GENOME ALIGNMENT INDEXES
        exon_db.__doc__ = 'Exon Annotation Database for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.dm2.exons', exon_db)
        msa.__doc__ = 'NLMSA Exon for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.NLMSA.dm2.exons', msa)
        exon_schema = pygr.Data.ManyToManyRelation(dm2, exon_db, bindAttrs = ('exon1',))
        exon_schema.__doc__ = 'Exon Schema for dm2'
        pygr.Data.addSchema('TEST.Annotation.NLMSA.dm2.exons', exon_schema)
        # BUILD ANNOTATION DATABASE FOR REFSEQ SPLICES
        splice_slices = pygr.Data.Collection(filename = os.path.join(self.path, 'refGene_spliceAnnot_dm2.cdb'), \
            intKeys = True, mode = 'c', writeback = False) # ONLY C
        splice_db = seqdb.AnnotationDB(splice_slices, dm2,
                               sliceAttrDict = dict(id = 0, splice_id = 1, orientation = 2,
                                                  gene_id = 3, start = 4, stop = 5))
        msa = cnestedlist.NLMSA(os.path.join(self.path, 'refGene_spliceAnnot_dm2'), 'w', \
            pairwiseMode = True, bidirectional = False)
        #for lines in open('refGene_spliceAnnot_dm2.txt', 'r').xreadlines():
        for lines in open('refGene_spliceAnnot_chrYh_dm2.txt', 'r').xreadlines():
            row = [x for x in lines.split('\t')] # CONVERT TO LIST SO MUTABLE
            row[1] = int(row[1]) # CONVERT FROM STRING TO INTEGER
            splice_slices[row[1]] = row
            splice = splice_db[row[1]] # GET THE ANNOTATION OBJECT FOR THIS EXON
            msa.addAnnotation(splice) # SAVE IT TO GENOME MAPPING
        splice_db.clear_cache() # not really necessary; cache should autoGC
        splice_slices.close() # SHELVE SHOULD BE EXPLICITLY CLOSED IN ORDER TO SAVE CURRENT CONTENTS
        msa.build() # FINALIZE GENOME ALIGNMENT INDEXES
        splice_db.__doc__ = 'Splice Annotation Database for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.dm2.splices', splice_db)
        msa.__doc__ = 'NLMSA Splice for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.NLMSA.dm2.splices', msa)
        splice_schema = pygr.Data.ManyToManyRelation(dm2, splice_db, bindAttrs = ('splice1',))
        splice_schema.__doc__ = 'Splice Schema for dm2'
        pygr.Data.addSchema('TEST.Annotation.NLMSA.dm2.splices', splice_schema)
        # BUILD ANNOTATION DATABASE FOR MOST CONSERVED ELEMENTS FROM UCSC
        ucsc_slices = pygr.Data.Collection(filename = os.path.join(self.path, 'phastConsElements15way_dm2.cdb'), \
            intKeys = True, mode = 'c', writeback = False) # ONLY C
        ucsc_db = seqdb.AnnotationDB(ucsc_slices, dm2,
                               sliceAttrDict = dict(id = 0, ucsc_id = 1, orientation = 2,
                                                  gene_id = 3, start = 4, stop = 5))
        msa = cnestedlist.NLMSA(os.path.join(self.path, 'phastConsElements15way_dm2'), 'w', \
            pairwiseMode = True, bidirectional = False)
        #for lines in open('phastConsElements15way_dm2.txt', 'r').xreadlines():
        for lines in open('phastConsElements15way_chrYh_dm2.txt', 'r').xreadlines():
            row = [x for x in lines.split('\t')] # CONVERT TO LIST SO MUTABLE
            row[1] = int(row[1]) # CONVERT FROM STRING TO INTEGER
            ucsc_slices[row[1]] = row
            ucsc = ucsc_db[row[1]] # GET THE ANNOTATION OBJECT FOR THIS EXON
            msa.addAnnotation(ucsc) # SAVE IT TO GENOME MAPPING
        ucsc_db.clear_cache() # not really necessary; cache should autoGC
        ucsc_slices.close() # SHELVE SHOULD BE EXPLICITLY CLOSED IN ORDER TO SAVE CURRENT CONTENTS
        msa.build() # FINALIZE GENOME ALIGNMENT INDEXES
        ucsc_db.__doc__ = 'Most Conserved Elements for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.UCSC.dm2.mostconserved', ucsc_db)
        msa.__doc__ = 'NLMSA for Most Conserved Elements for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.UCSC.NLMSA.dm2.mostconserved', msa)
        ucsc_schema = pygr.Data.ManyToManyRelation(dm2, ucsc_db, bindAttrs = ('element1',))
        ucsc_schema.__doc__ = 'Schema for UCSC Most Conserved Elements for dm2'
        pygr.Data.addSchema('TEST.Annotation.UCSC.NLMSA.dm2.mostconserved', ucsc_schema)
        pygr.Data.save()
        reload(pygr.Data)

        # QUERY TO EXON AND SPLICES ANNOTATION DATABASE
        dm2 = pygr.Data.getResource('TEST.Seq.Genome.dm2')
        exonmsa = pygr.Data.getResource('TEST.Annotation.NLMSA.dm2.exons')
        splicemsa = pygr.Data.getResource('TEST.Annotation.NLMSA.dm2.splices')
        conservedmsa = pygr.Data.getResource('TEST.Annotation.UCSC.NLMSA.dm2.mostconserved')
        exons = pygr.Data.getResource('TEST.Annotation.dm2.exons')
        splices = pygr.Data.getResource('TEST.Annotation.dm2.splices')
        mostconserved = pygr.Data.getResource('TEST.Annotation.UCSC.dm2.mostconserved')

        # OPEN DM2_MULTIZ15WAY NLMSA
        msa = cnestedlist.NLMSA(os.path.join(msaDir, 'dm2_multiz15way'), 'r', trypath = [seqDir])

        #exonAnnotFileName = 'Annotation_ConservedElement_Exons_dm2.txt'
        #intronAnnotFileName = 'Annotation_ConservedElement_Introns_dm2.txt'
        exonAnnotFileName = 'Annotation_ConservedElement_Exons_chrYh_dm2.txt' # FOR TESTING
        intronAnnotFileName = 'Annotation_ConservedElement_Introns_chrYh_dm2.txt' # FOR TESTING
        newexonAnnotFileName = os.path.join(self.path, 'new_Exons_dm2.txt')
        newintronAnnotFileName = os.path.join(self.path, 'new_Introns_dm2.txt')
        tmpexonAnnotFileName = self.copyFile(exonAnnotFileName)
        tmpintronAnnotFileName = self.copyFile(intronAnnotFileName)

        chrList = dm2.seqLenDict.keys()
        chrList.sort()
        chrList = ['chrYh'] # FOR TESTING

        outfile = open(newexonAnnotFileName, 'w')
        for chrid in chrList:
            slice = dm2[chrid]
            try:
                ex1 = exonmsa[slice]
            except KeyError:
                continue
            else:
                exlist1 = [(ix.exon_id, ix) for ix in ex1.keys()]
                exlist1.sort()
                for ixx, exon in exlist1:
                    saveList = []
                    tmp = exon.sequence
                    tmpexon = exons[exon.exon_id]
                    tmpslice = tmpexon.sequence # FOR REAL EXON COORDINATE
                    wlist1 = 'EXON', chrid, tmpexon.exon_id, tmpexon.gene_id, tmpslice.start, tmpslice.stop
                    try:
                        out1 = conservedmsa[tmp]
                    except KeyError:
                        pass
                    else:
                        elementlist = [(ix.ucsc_id, ix) for ix in out1.keys()]
                        elementlist.sort()
                        for iyy, element in elementlist:
                            if element.stop - element.start < 100: continue
                            score = int(string.split(element.gene_id, '=')[1])
                            if score < 100: continue
                            tmp2 = element.sequence
                            tmpelement = mostconserved[element.ucsc_id]
                            tmpslice2 = tmpelement.sequence # FOR REAL ELEMENT COORDINATE
                            wlist2 = wlist1 + (tmpelement.ucsc_id, tmpelement.gene_id, tmpslice2.start, tmpslice2.stop)
                            slicestart, sliceend = max(tmp.start, tmp2.start), min(tmp.stop, tmp2.stop)
                            tmp1 = msa.seqDict['dm2.' + chrid][slicestart:sliceend]
                            edges = msa[tmp1].edges()
                            for src, dest, e in edges:
                                if src.stop - src.start < 100: continue
                                palign, pident = e.pAligned(), e.pIdentity()
                                if palign < 0.8 or pident < 0.8: continue
                                palign, pident = '%.2f' % palign, '%.2f' % pident
                                wlist3 = wlist2 + ((~msa.seqDict)[src], str(src), src.start, src.stop, \
                                    (~msa.seqDict)[dest], \
                                    str(dest), dest.start, dest.stop, palign, pident)
                                saveList.append('\t'.join(map(str, wlist3)) + '\n')
                        saveList.sort()
                        for saveline in saveList:
                            outfile.write(saveline)
        outfile.close()
        import md5
        md5old = md5.new()
        md5old.update(open(tmpexonAnnotFileName, 'r').read())
        md5new = md5.new()
        md5new.update(open(newexonAnnotFileName, 'r').read())
        assert md5old.digest() == md5new.digest() # MD5 COMPARISON INSTEAD OF COMPARING EACH CONTENTS

        outfile = open(newintronAnnotFileName, 'w')
        for chrid in chrList:
            slice = dm2[chrid]
            try:
                sp1 = splicemsa[slice]
            except:
                continue
            else:
                splist1 = [(ix.splice_id, ix) for ix in sp1.keys()]
                splist1.sort()
                for ixx, splice in splist1:
                    saveList = []
                    tmp = splice.sequence
                    tmpsplice = splices[splice.splice_id]
                    tmpslice = tmpsplice.sequence # FOR REAL EXON COORDINATE
                    wlist1 = 'INTRON', chrid, tmpsplice.splice_id, tmpsplice.gene_id, tmpslice.start, tmpslice.stop
                    try:
                        out1 = conservedmsa[tmp]
                    except KeyError:
                        pass
                    else:
                        elementlist = [(ix.ucsc_id, ix) for ix in out1.keys()]
                        elementlist.sort()
                        for iyy, element in elementlist:
                            if element.stop - element.start < 100: continue
                            score = int(string.split(element.gene_id, '=')[1])
                            if score < 100: continue
                            tmp2 = element.sequence
                            tmpelement = mostconserved[element.ucsc_id]
                            tmpslice2 = tmpelement.sequence # FOR REAL ELEMENT COORDINATE
                            wlist2 = wlist1 + (tmpelement.ucsc_id, tmpelement.gene_id, tmpslice2.start, tmpslice2.stop)
                            slicestart, sliceend = max(tmp.start, tmp2.start), min(tmp.stop, tmp2.stop)
                            tmp1 = msa.seqDict['dm2.' + chrid][slicestart:sliceend]
                            edges = msa[tmp1].edges()
                            for src, dest, e in edges:
                                if src.stop - src.start < 100: continue
                                palign, pident = e.pAligned(), e.pIdentity()
                                if palign < 0.8 or pident < 0.8: continue
                                palign, pident = '%.2f' % palign, '%.2f' % pident
                                wlist3 = wlist2 + ((~msa.seqDict)[src], str(src), src.start, src.stop, \
                                    (~msa.seqDict)[dest], \
                                    str(dest), dest.start, dest.stop, palign, pident)
                                saveList.append('\t'.join(map(str, wlist3)) + '\n')
                        saveList.sort()
                        for saveline in saveList:
                            outfile.write(saveline)
        outfile.close()
        import md5
        md5old = md5.new()
        md5old.update(open(tmpintronAnnotFileName, 'r').read())
        md5new = md5.new()
        md5new.update(open(newintronAnnotFileName, 'r').read())
        assert md5old.digest() == md5new.digest() # MD5 COMPARISON INSTEAD OF COMPARING EACH CONTENTS

    def mysqlannot_test(self): # BUILD ANNOTATION DB FROM MYSQL
        os.environ['PYGRDATAPATH'] = self.path
        import pygr.Data
        from pygr import seqdb, cnestedlist, sqlgraph
        dm2 = pygr.Data.getResource('TEST.Seq.Genome.dm2')
        # BUILD ANNOTATION DATABASE FOR REFSEQ EXONS: MYSQL VERSION
        #exon_slices = sqlgraph.SQLTableClustered('PYGRDB_JAN06.pygr_refGene_exonAnnot_dm2',
        exon_slices = sqlgraph.SQLTableClustered('PYGRDB_JAN06.pygr_refGene_exonAnnot_chrYh_dm2',
            clusterKey = 'chromosome', maxCache = 0)
        exon_db = seqdb.AnnotationDB(exon_slices, dm2, sliceAttrDict = dict(id = 'chromosome', \
            gene_id = 'name', exon_id = 'exon_id'))
        msa = cnestedlist.NLMSA(os.path.join(self.path, 'refGene_exonAnnot_SQL_dm2'), 'w', \
            pairwiseMode = True, bidirectional = False)
        for id in exon_db:
            msa.addAnnotation(exon_db[id])
        exon_db.clear_cache() # not really necessary; cache should autoGC
        exon_slices.clear_cache()
        msa.build()
        exon_db.__doc__ = 'SQL Exon Annotation Database for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.SQL.dm2.exons', exon_db)
        msa.__doc__ = 'SQL NLMSA Exon for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.NLMSA.SQL.dm2.exons', msa)
        exon_schema = pygr.Data.ManyToManyRelation(dm2, exon_db, bindAttrs = ('exon2',))
        exon_schema.__doc__ = 'SQL Exon Schema for dm2'
        pygr.Data.addSchema('TEST.Annotation.NLMSA.SQL.dm2.exons', exon_schema)
        # BUILD ANNOTATION DATABASE FOR REFSEQ SPLICES: MYSQL VERSION
        #splice_slices = sqlgraph.SQLTableClustered('PYGRDB_JAN06.pygr_refGene_spliceAnnot_dm2',
        splice_slices = sqlgraph.SQLTableClustered('PYGRDB_JAN06.pygr_refGene_spliceAnnot_chrYh_dm2',
            clusterKey = 'chromosome', maxCache = 0)
        splice_db = seqdb.AnnotationDB(splice_slices, dm2, sliceAttrDict = dict(id = 'chromosome', \
            gene_id = 'name', splice_id = 'splice_id'))
        msa = cnestedlist.NLMSA(os.path.join(self.path, 'refGene_spliceAnnot_SQL_dm2'), 'w', \
            pairwiseMode = True, bidirectional = False)
        for id in splice_db:
            msa.addAnnotation(splice_db[id])
        splice_db.clear_cache() # not really necessary; cache should autoGC
        splice_slices.clear_cache()
        msa.build()
        splice_db.__doc__ = 'SQL Splice Annotation Database for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.SQL.dm2.splices', splice_db)
        msa.__doc__ = 'SQL NLMSA Splice for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.NLMSA.SQL.dm2.splices', msa)
        splice_schema = pygr.Data.ManyToManyRelation(dm2, splice_db, bindAttrs = ('splice2',))
        splice_schema.__doc__ = 'SQL Splice Schema for dm2'
        pygr.Data.addSchema('TEST.Annotation.NLMSA.SQL.dm2.splices', splice_schema)
        # BUILD ANNOTATION DATABASE FOR MOST CONSERVED ELEMENTS FROM UCSC: MYSQL VERSION
        #ucsc_slices = sqlgraph.SQLTableClustered('PYGRDB_JAN06.pygr_phastConsElements15way_dm2',
        ucsc_slices = sqlgraph.SQLTableClustered('PYGRDB_JAN06.pygr_phastConsElements15way_chrYh_dm2',
            clusterKey = 'chromosome', maxCache = 0)
        ucsc_db = seqdb.AnnotationDB(ucsc_slices, dm2, sliceAttrDict = dict(id = 'chromosome', \
            gene_id = 'name', ucsc_id = 'ucsc_id'))
        msa = cnestedlist.NLMSA(os.path.join(self.path, 'phastConsElements15way_SQL_dm2'), 'w', \
            pairwiseMode = True, bidirectional = False)
        for id in ucsc_db:
            msa.addAnnotation(ucsc_db[id])
        ucsc_db.clear_cache() # not really necessary; cache should autoGC
        ucsc_slices.clear_cache()
        msa.build()
        ucsc_db.__doc__ = 'SQL Most Conserved Elements for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.UCSC.SQL.dm2.mostconserved', ucsc_db)
        msa.__doc__ = 'SQL NLMSA for Most Conserved Elements for dm2'
        pygr.Data.getResource.addResource('TEST.Annotation.UCSC.NLMSA.SQL.dm2.mostconserved', msa)
        ucsc_schema = pygr.Data.ManyToManyRelation(dm2, ucsc_db, bindAttrs = ('element2',))
        ucsc_schema.__doc__ = 'SQL Schema for UCSC Most Conserved Elements for dm2'
        pygr.Data.addSchema('TEST.Annotation.UCSC.NLMSA.SQL.dm2.mostconserved', ucsc_schema)
        pygr.Data.save()
        reload(pygr.Data)

        # QUERY TO EXON AND SPLICES ANNOTATION DATABASE
        dm2 = pygr.Data.getResource('TEST.Seq.Genome.dm2')
        exonmsa = pygr.Data.getResource('TEST.Annotation.NLMSA.SQL.dm2.exons')
        splicemsa = pygr.Data.getResource('TEST.Annotation.NLMSA.SQL.dm2.splices')
        conservedmsa = pygr.Data.getResource('TEST.Annotation.UCSC.NLMSA.SQL.dm2.mostconserved')
        exons = pygr.Data.getResource('TEST.Annotation.SQL.dm2.exons')
        splices = pygr.Data.getResource('TEST.Annotation.SQL.dm2.splices')
        mostconserved = pygr.Data.getResource('TEST.Annotation.UCSC.SQL.dm2.mostconserved')

        # OPEN DM2_MULTIZ15WAY NLMSA
        msa = cnestedlist.NLMSA(os.path.join(msaDir, 'dm2_multiz15way'), 'r', trypath = [seqDir])

        #exonAnnotFileName = 'Annotation_ConservedElement_Exons_dm2.txt'
        #intronAnnotFileName = 'Annotation_ConservedElement_Introns_dm2.txt'
        exonAnnotFileName = 'Annotation_ConservedElement_Exons_chrYh_dm2.txt' # FOR TESTING
        intronAnnotFileName = 'Annotation_ConservedElement_Introns_chrYh_dm2.txt' # FOR TESTING
        newexonAnnotFileName = os.path.join(self.path, 'new_Exons_dm2.txt')
        newintronAnnotFileName = os.path.join(self.path, 'new_Introns_dm2.txt')
        tmpexonAnnotFileName = self.copyFile(exonAnnotFileName)
        tmpintronAnnotFileName = self.copyFile(intronAnnotFileName)

        chrList = dm2.seqLenDict.keys()
        chrList.sort()
        chrList = ['chrYh'] # FOR TESTING

        outfile = open(newexonAnnotFileName, 'w')
        for chrid in chrList:
            slice = dm2[chrid]
            try:
                ex1 = exonmsa[slice]
            except KeyError:
                continue
            else:
                exlist1 = [(ix.exon_id, ix) for ix in ex1.keys()]
                exlist1.sort()
                for ixx, exon in exlist1:
                    saveList = []
                    tmp = exon.sequence
                    tmpexon = exons[exon.exon_id]
                    tmpslice = tmpexon.sequence # FOR REAL EXON COORDINATE
                    wlist1 = 'EXON', chrid, tmpexon.exon_id, tmpexon.gene_id, tmpslice.start, tmpslice.stop
                    try:
                        out1 = conservedmsa[tmp]
                    except KeyError:
                        pass
                    else:
                        elementlist = [(ix.ucsc_id, ix) for ix in out1.keys()]
                        elementlist.sort()
                        for iyy, element in elementlist:
                            if element.stop - element.start < 100: continue
                            score = int(string.split(element.gene_id, '=')[1])
                            if score < 100: continue
                            tmp2 = element.sequence
                            tmpelement = mostconserved[element.ucsc_id]
                            tmpslice2 = tmpelement.sequence # FOR REAL ELEMENT COORDINATE
                            wlist2 = wlist1 + (tmpelement.ucsc_id, tmpelement.gene_id, tmpslice2.start, tmpslice2.stop)
                            slicestart, sliceend = max(tmp.start, tmp2.start), min(tmp.stop, tmp2.stop)
                            tmp1 = msa.seqDict['dm2.' + chrid][slicestart:sliceend]
                            edges = msa[tmp1].edges()
                            for src, dest, e in edges:
                                if src.stop - src.start < 100: continue
                                palign, pident = e.pAligned(), e.pIdentity()
                                if palign < 0.8 or pident < 0.8: continue
                                palign, pident = '%.2f' % palign, '%.2f' % pident
                                wlist3 = wlist2 + ((~msa.seqDict)[src], str(src), src.start, src.stop, \
                                    (~msa.seqDict)[dest], \
                                    str(dest), dest.start, dest.stop, palign, pident)
                                saveList.append('\t'.join(map(str, wlist3)) + '\n')
                        saveList.sort()
                        for saveline in saveList:
                            outfile.write(saveline)
        outfile.close()
        import md5
        md5old = md5.new()
        md5old.update(open(tmpexonAnnotFileName, 'r').read())
        md5new = md5.new()
        md5new.update(open(newexonAnnotFileName, 'r').read())
        assert md5old.digest() == md5new.digest() # MD5 COMPARISON INSTEAD OF COMPARING EACH CONTENTS

        outfile = open(newintronAnnotFileName, 'w')
        for chrid in chrList:
            slice = dm2[chrid]
            try:
                sp1 = splicemsa[slice]
            except:
                continue
            else:
                splist1 = [(ix.splice_id, ix) for ix in sp1.keys()]
                splist1.sort()
                for ixx, splice in splist1:
                    saveList = []
                    tmp = splice.sequence
                    tmpsplice = splices[splice.splice_id]
                    tmpslice = tmpsplice.sequence # FOR REAL EXON COORDINATE
                    wlist1 = 'INTRON', chrid, tmpsplice.splice_id, tmpsplice.gene_id, tmpslice.start, tmpslice.stop
                    try:
                        out1 = conservedmsa[tmp]
                    except KeyError:
                        pass
                    else:
                        elementlist = [(ix.ucsc_id, ix) for ix in out1.keys()]
                        elementlist.sort()
                        for iyy, element in elementlist:
                            if element.stop - element.start < 100: continue
                            score = int(string.split(element.gene_id, '=')[1])
                            if score < 100: continue
                            tmp2 = element.sequence
                            tmpelement = mostconserved[element.ucsc_id]
                            tmpslice2 = tmpelement.sequence # FOR REAL ELEMENT COORDINATE
                            wlist2 = wlist1 + (tmpelement.ucsc_id, tmpelement.gene_id, tmpslice2.start, tmpslice2.stop)
                            slicestart, sliceend = max(tmp.start, tmp2.start), min(tmp.stop, tmp2.stop)
                            tmp1 = msa.seqDict['dm2.' + chrid][slicestart:sliceend]
                            edges = msa[tmp1].edges()
                            for src, dest, e in edges:
                                if src.stop - src.start < 100: continue
                                palign, pident = e.pAligned(), e.pIdentity()
                                if palign < 0.8 or pident < 0.8: continue
                                palign, pident = '%.2f' % palign, '%.2f' % pident
                                wlist3 = wlist2 + ((~msa.seqDict)[src], str(src), src.start, src.stop, \
                                    (~msa.seqDict)[dest], \
                                    str(dest), dest.start, dest.stop, palign, pident)
                                saveList.append('\t'.join(map(str, wlist3)) + '\n')
                        saveList.sort()
                        for saveline in saveList:
                            outfile.write(saveline)
        outfile.close()
        import md5
        md5old = md5.new()
        md5old.update(open(tmpintronAnnotFileName, 'r').read())
        md5new = md5.new()
        md5new.update(open(newintronAnnotFileName, 'r').read())
        assert md5old.digest() == md5new.digest() # MD5 COMPARISON INSTEAD OF COMPARING EACH CONTENTS

