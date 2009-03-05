
import sys, os, string, glob

mafDir = '/result/pygr_megatest/maf_data'
seqDir = '/result/pygr_megatest/seq_data'

## mafDir CONTAINS FOLLOWING DM2 MULTIZ15WAY MAF ALIGNMENTS
## seqDir CONTAINS FOLLOWING 15 GENOME ASSEMBLIES AND THEIR SEQDB FILES
## TEST INPUT/OUPTUT FOR COMPARISON, THESE FILES SHOULD BE IN THIS DIRECTORY
##        outfileName = 'splicesite_dm2.txt' # CHR4H TESTING
##        outputName = 'splicesite_dm2_multiz15way.txt' # CHR4H TESTING
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
    def build_test(self): # BUILD NLMSA AND QUERY RESULT COMPARISON
        os.environ['PYGRDATAPATH'] = self.path
        import pygr.Data
        from pygr import seqdb, cnestedlist
        genomedict = {}
        for orgstr in msaSpeciesList:
            genomedict[orgstr] = pygr.Data.getResource('TEST.Seq.Genome.' + orgstr)
        uniondict = seqdb.PrefixUnionDict(genomedict)
        import glob
        maflist = glob.glob(os.path.join(mafDir, 'chr4h.maf')) # CHR4H TESTING
        maflist.sort()
        msaname = os.path.join(self.path, 'dm2_multiz15way')
        msa1 = cnestedlist.NLMSA(msaname, 'w', uniondict, maflist, maxlen = 536870912, maxint = 22369620) # 500MB VERSION
        msa1.save_seq_dict()
        msa1.__doc__ = 'TEST NLMSA for dm2 multiz15way'
        pygr.Data.getResource.addResource('TEST.MSA.UCSC.dm2_multiz15way', msa1)
        pygr.Data.save()
        msa = pygr.Data.getResource('TEST.MSA.UCSC.dm2_multiz15way')
        outfileName = 'splicesite_dm2_chr4h.txt' # CHR4H TESTING
        outputName = 'splicesite_dm2_chr4h_multiz15way.txt' # CHR4H TESTING
        newOutputName = 'splicesite_new1.txt'
        tmpInputName = self.copyFile(outfileName)
        tmpOutputName = self.copyFile(outputName)
        tmpNewOutputName = os.path.join(self.path, newOutputName)
        outfile = open(tmpNewOutputName, 'w')
        for lines in open(tmpInputName, 'r').xreadlines():
            chrid, intstart, intend, nobs = string.split(lines.strip(), '\t')
            intstart, intend, nobs = int(intstart), int(intend), int(nobs)
            site1 = msa.seqDict['dm2' + '.' + chrid][intstart:intstart+2]
            site2 = msa.seqDict['dm2' + '.' + chrid][intend-2:intend]
            edges1 = msa[site1].edges()
            edges2 = msa[site2].edges()
            if len(edges1) == 0: # EMPTY EDGES
                wlist = str(site1), 'dm2', chrid, intstart, intstart+2, '', '', '', '', ''
                outfile.write('\t'.join(map(str, wlist)) + '\n')
            if len(edges2) == 0: # EMPTY EDGES
                wlist = str(site2), 'dm2', chrid, intend-2, intend, '', '', '', '', ''
                outfile.write('\t'.join(map(str, wlist)) + '\n')
            saveList = []
            for src, dest, e in edges1:
                if len(str(src)) != 2 or len(str(dest)) != 2: continue
                dotindex = (~msa.seqDict)[src].index('.')
                srcspecies, src1 = (~msa.seqDict)[src][:dotindex], (~msa.seqDict)[src][dotindex+1:]
                dotindex = (~msa.seqDict)[dest].index('.')
                destspecies, dest1 = (~msa.seqDict)[dest][:dotindex], (~msa.seqDict)[dest][dotindex+1:]
                wlist = str(src), srcspecies, src1, src.start, src.stop, str(dest), \
                    destspecies, dest1, dest.start, dest.stop
                saveList.append('\t'.join(map(str, wlist)) + '\n')
            for src, dest, e in edges2:
                if len(str(src)) != 2 or len(str(dest)) != 2: continue
                dotindex = (~msa.seqDict)[src].index('.')
                srcspecies, src1 = (~msa.seqDict)[src][:dotindex], (~msa.seqDict)[src][dotindex+1:]
                dotindex = (~msa.seqDict)[dest].index('.')
                destspecies, dest1 = (~msa.seqDict)[dest][:dotindex], (~msa.seqDict)[dest][dotindex+1:]
                wlist = str(src), srcspecies, src1, src.start, src.stop, str(dest), \
                    destspecies, dest1, dest.start, dest.stop
                saveList.append('\t'.join(map(str, wlist)) + '\n')
            saveList.sort() # SORTED IN ORDER TO COMPARE WITH PREVIOUS RESULTS
            for saveline in saveList:
                outfile.write(saveline)
        outfile.close()
        import md5
        md5old = md5.new()
        md5old.update(open(tmpNewOutputName, 'r').read())
        md5new = md5.new()
        md5new.update(open(tmpOutputName, 'r').read())
        assert md5old.digest() == md5new.digest() # MD5 COMPARISON INSTEAD OF COMPARING EACH CONTENTS

        # TEXT<->BINARY TEST
        msafilelist = glob.glob(msaname + '*')
        cnestedlist.dump_textfile(msaname, os.path.join(self.path, 'dm2_multiz15way.txt'))
        for filename in msafilelist: os.remove(filename)
        runPath = os.path.realpath(os.curdir)
        os.chdir(self.path)
        cnestedlist.textfile_to_binaries('dm2_multiz15way.txt')
        os.chdir(runPath)

        msa1 = cnestedlist.NLMSA(msaname, 'r')
        msa1.__doc__ = 'TEST NLMSA for dm2 multiz15way'
        pygr.Data.getResource.addResource('TEST.MSA.UCSC.dm2_multiz15way', msa1)
        pygr.Data.save()
        msa = pygr.Data.getResource('TEST.MSA.UCSC.dm2_multiz15way')
        newOutputName = 'splicesite_new2.txt'
        tmpInputName = self.copyFile(outfileName)
        tmpOutputName = self.copyFile(outputName)
        tmpNewOutputName = os.path.join(self.path, newOutputName)
        outfile = open(tmpNewOutputName, 'w')
        for lines in open(tmpInputName, 'r').xreadlines():
            chrid, intstart, intend, nobs = string.split(lines.strip(), '\t')
            intstart, intend, nobs = int(intstart), int(intend), int(nobs)
            site1 = msa.seqDict['dm2' + '.' + chrid][intstart:intstart+2]
            site2 = msa.seqDict['dm2' + '.' + chrid][intend-2:intend]
            edges1 = msa[site1].edges()
            edges2 = msa[site2].edges()
            if len(edges1) == 0: # EMPTY EDGES
                wlist = str(site1), 'dm2', chrid, intstart, intstart+2, '', '', '', '', ''
                outfile.write('\t'.join(map(str, wlist)) + '\n')
            if len(edges2) == 0: # EMPTY EDGES
                wlist = str(site2), 'dm2', chrid, intend-2, intend, '', '', '', '', ''
                outfile.write('\t'.join(map(str, wlist)) + '\n')
            saveList = []
            for src, dest, e in edges1:
                if len(str(src)) != 2 or len(str(dest)) != 2: continue
                dotindex = (~msa.seqDict)[src].index('.')
                srcspecies, src1 = (~msa.seqDict)[src][:dotindex], (~msa.seqDict)[src][dotindex+1:]
                dotindex = (~msa.seqDict)[dest].index('.')
                destspecies, dest1 = (~msa.seqDict)[dest][:dotindex], (~msa.seqDict)[dest][dotindex+1:]
                wlist = str(src), srcspecies, src1, src.start, src.stop, str(dest), \
                    destspecies, dest1, dest.start, dest.stop
                saveList.append('\t'.join(map(str, wlist)) + '\n')
            for src, dest, e in edges2:
                if len(str(src)) != 2 or len(str(dest)) != 2: continue
                dotindex = (~msa.seqDict)[src].index('.')
                srcspecies, src1 = (~msa.seqDict)[src][:dotindex], (~msa.seqDict)[src][dotindex+1:]
                dotindex = (~msa.seqDict)[dest].index('.')
                destspecies, dest1 = (~msa.seqDict)[dest][:dotindex], (~msa.seqDict)[dest][dotindex+1:]
                wlist = str(src), srcspecies, src1, src.start, src.stop, str(dest), \
                    destspecies, dest1, dest.start, dest.stop
                saveList.append('\t'.join(map(str, wlist)) + '\n')
            saveList.sort() # SORTED IN ORDER TO COMPARE WITH PREVIOUS RESULTS
            for saveline in saveList:
                outfile.write(saveline)
        outfile.close()
        import md5
        md5old = md5.new()
        md5old.update(open(tmpNewOutputName, 'r').read())
        md5new = md5.new()
        md5new.update(open(tmpOutputName, 'r').read())
        assert md5old.digest() == md5new.digest() # MD5 COMPARISON INSTEAD OF COMPARING EACH CONTENTS

