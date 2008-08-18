''' # These two classes are replaced by the magic sqlgraph.SQLGraph class
class ExonToTranscriptInv(object):
    'inverse of a TranscriptToExon mapping: Ensembl exon obj -> Ensembl transcript obj'

    def __init__(self, mapper):
        self.inverseDB = mapper
    
    def __getitem__(self, k):
        'find corresponding transcripts to the given exon'

        exonID = k.id
        n = self.inverseDB.cursor.execute('select transcript_id from %s where exon_id = %s' %(self.inverseDB.exon_transcript, exonID))
        t = self.inverseDB.cursor.fetchall()
        transcripts = []
        for row in t:
            id = row[0]
            transcript = self.inverseDB.transcriptDB[id]
            transcripts.append(transcript)
        return transcripts

    def __invert__(self):
        return self.inverseDB

class ExonToTranscript(object):
    'Provide a mapping of a transcript obj -> a set of exon objects'

    def __init__(self, transcriptDB, exonDB):
        self.exonDB = exonDB
        self.transcriptDB = transcriptDB
        self.inverseDB = TranscriptToExonInv(self)
        self.exon_transcript = self.exonDB.name.split('.')[0] + '.exon_transcript'
        self.cursor = self.exonDB.cursor 

    def __getitem__(self, k):
        transcriptID = k.id
        n = self.cursor.execute('select exon_id from %s where transcript_id = %s' %(self.exon_transcript, transcriptID))
        t = self.cursor.fetchall()
        exons = []
        for row in t:
            id = row[0]
            exon = self.exonDB[id]
            exons.append(exon)
        return exons

    def __invert__(self):
        return self.inverseDB
'''
class GeneToTranscriptInv(object):
    'inverse of a GeneToTranscript mapping: an Ensembl transcript obj -> an Ensembl gene obj'

    def __init__(self, mapper):
        self.inverseDB = mapper
    
    def __getitem__(self, k):
        'find the corresponding gene to the given transcript'

        geneID = k.gene_id
        gene = self.inverseDB.geneDB[geneID] 
        return gene

    def __invert__(self):
        return self.inverseDB

class GeneToTranscript(object):
    'Provide a mapping of a gene obj -> a set of transcript objects'

    def __init__(self, geneDB, transcriptDB):
        self.geneDB = geneDB
        self.transcriptDB = transcriptDB
        self.inverseDB = GeneToTranscriptInv(self)

    def __getitem__(self, k):
        geneID = k.id
        #n = self.cursor.execute('select exon_id from %s where transcript_id = %s' %(self.exon_transcript, transcriptID))
        t = self.transcriptDB.select('where gene_id = %s', (geneID))
        transcripts = []
        for row in t:
            #id = row[0]
            #exon = self.exonDB[id]
            transcripts.append(row)
        return transcripts

    def __invert__(self):
        return self.inverseDB


class TranscriptToTranslationInv(object):
    'inverse of a TranscriptToTranslation mapping: an Ensembl translation obj -> an Ensembl transcript obj'

    def __init__(self, mapper):
        self.inverseDB = mapper
    
    def __getitem__(self, k):
        'find the corresponding transcript to the given translation'

        transcriptID = k.transcript_id
        transcript = self.inverseDB.transcriptDB[transcriptID] 
        return transcript

    def __invert__(self):
        return self.inverseDB

class TranscriptToTranslation(object):
    'Provide a mapping of a transcript obj -> a translation object'

    def __init__(self, transcriptDB, translationDB):
        self.transcriptDB = transcriptDB
        self.translationDB = translationDB
        self.inverseDB = TranscriptToTranslationInv(self)

    def __getitem__(self, k):
        transcriptID = k.id
        t = self.translationDB.select('where transcript_id = %s', (transcriptID))
        translation = None
        translations = []
        for row in t:
            #id = row[0]
            #exon = self.exonDB[id]
            translations.append(row)
        n = len(translations)
        if n>1:
            raise KeyError('duplicated! Transcript -> Translation is not one to one mapping!')
        if n == 0:
            print 'This transcript is not translateable!'
        else:
            translation = translations[0]
        return translation

    def __invert__(self):
        return self.inverseDB





class PtranscriptToPexonInv(object):
    'inverse of a PtranscriptToPexon mapping: an Ensembl prediction_exon obj -> an Ensembl prediction_transcript obj'

    def __init__(self, mapper):
        self.inverseDB = mapper
    
    def __getitem__(self, k):
        'find the corresponding transcript to the given exon'

        ptranscriptID = k.prediction_transcript_id
        ptranscript = self.inverseDB.ptranscriptDB[ptranscriptID] 
        return ptranscript

    def __invert__(self):
        return self.inverseDB

class PtranscriptToPexon(object):
    'Provide a mapping of a prediction_transcript obj -> a set of prediction_exon objects'

    def __init__(self, ptranscriptDB, pexonDB):
        self.pexonDB = pexonDB
        self.ptranscriptDB = ptranscriptDB
        self.inverseDB = PtranscriptToPexonInv(self)
        #self.cursor = self.ptranscriptDB.cursor 

    def __getitem__(self, k):
        ptranscriptID = k.id
        #n = self.cursor.execute('select exon_id from %s where transcript_id = %s' %(self.exon_transcript, transcriptID))
        t = self.pexonDB.select('where prediction_transcript_id = %s', (ptranscriptID))
        #t = self.cursor.fetchall()
        pexons = []
        for row in t:
            #id = row[0]
            #exon = self.exonDB[id]
            pexons.append(row)
        return pexons

    def __invert__(self):
        return self.inverseDB
