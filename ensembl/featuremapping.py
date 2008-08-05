class TranscriptToExonInv(object):
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

class TranscriptToExon(object):
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
