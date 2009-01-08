import pygr.Data
from ensembl import adaptor
from pygr import sqlgraph

def sample_DB():
    """sample code on how to retrieve ensembl features and schemas using a perl-like interface""" 

    serverRegistry = adaptor.get_registry(host='ensembldb.ensembl.org', user='anonymous')
    coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    # get feature annotation DB
    exonDB = coreDBAdaptor.get_feature('exon')
    exon = exonDB[1]
    print str(exon.sequence)
    transcriptDB = coreDBAdaptor.get_feature('transcript')
    transcript = transcriptDB[1]
    print 'length of this transcript: ', len(transcript.sequence)
    geneDB = coreDBAdaptor.get_feature('gene')
    gene = geneDB[1]
    
    # Retrieve annotation objects through their related annotation objects.  In the meantime, save the inter-relations into pygr.Data
    exons = transcript.get_all_exons()
    for e in exons:
        print e.id, e.seq_region_start, len(e.sequence)

    transcripts = exon.get_all_transcripts()
    for t in transcripts:
        print t.id, t.seq_region_start, len(t.sequence)

    genes = transcript.get_gene()
    for g in genes:
        print g.id, len(g.sequence)

    transcripts = gene.get_all_transcripts()
    for t in transcripts:
        print t.id, len(t.sequence)


    # get graph of feature annotation DBs
    #exonTranscript = coreDBAdaptor._create_graph(exonDB, transcriptDB, 'exon_transcript')
    #coreDBAdaptor._save_graph(exonTranscript, exonDB, transcriptDB, 'ManyToMany', 'transcripts', 'exons')

    #geneTranscript = coreDBAdaptor._create_graph(geneDB, transcriptDB, 'transcript')
    #coreDBAdaptor._save_graph(geneTranscript, geneDB, transcriptDB, 'OneToMany', 'transcripts', 'gene')
    #conn = coreDBAdaptor.conn
    #exonTranscript = sqlgraph.SQLGraph('homo_sapiens_core_47_36i.exon_transcript', serverInfo=conn, sourceDB=exonDB, targetDB=transcriptDB, attrAlias=dict(source_id='exon_id', target_id='transcript_id'))
    #exonTranscript.__doc__= 'an ensembl graph exonAnnoDB -> transcriptAnnoDB (homo_sapiens_core_47_36i)'
    #pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.exon_transcript = exonTranscript
    #pygr.Data.schema.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.exon_transcript = pygr.Data.ManyToManyRelation(exonDB, transcriptDB, bindAttrs=('transcripts', 'exons'))
    #pygr.Data.save()

def sample_pygrData():
    
    # get feature annotationDB from pygr.Data
    exonDB = pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.exon()
    transcriptDB = pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.transcript()
    geneDB = pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.gene()
    exon = exonDB[1]
    transcript = transcriptDB[1]
    gene = geneDB[1]

    # get feature annotation graphs from pygr.Data
    exonTranscript = pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.exon_transcript()
    geneTranscript = pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.gene_transcript()
    
    print '\n'
    exons = (~exonTranscript)[transcript]
    for e in exons:
        print e.id, e.seq_region_start, len(e.sequence)
    
    transcripts = exonTranscript[exon]
    for t in transcripts:
        print t.id, t.seq_region_start, len(t.sequence)

    genes = (~geneTranscript)[transcript]
    for g in genes:
        print g.id, len(g.sequence)

    transcripts = geneTranscript[gene]
    for t in transcripts:
        print t.id, len(t.sequence)

    #exonTranscript = pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.exon_transcript()
    #exons = (~exonTranscript)[transcript]
    #for e in exons:
    #    print e.id, len(e.sequence)
    #transcript = exonTranscript[exon]
    #for t in transcript:
    #    print t.id, len(t.sequence)
    #print len(gene.sequence)
    #transcript = transcriptDB[2]
    #transcript = transcriptDB[1]
    #geneTranscript = pygr.Data.Bio.Annotation.Ensembl.homo_sapiens_core_47_36i.gene_transcript()
    #transcripts = geneTranscript[gene]
    #for t in transcripts:
    #    print t.id, len(t.sequence)
    #gene = (~geneTranscript)[transcript]
    #for g in gene:
    #    print g.id, len(g.sequence)

if __name__ == '__main__': # example code
    
    sample_DB()
    sample_pygrData()

    '''
    serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
    coreDBAdaptor = serverRegistry.get_DBAdaptor('homo_sapiens', 'core', '47_36i')
    exonDB = coreDBAdaptor.get_feature('exon')
    exon = exonDB[1]
    print str(exon.sequence)
    #slice = coreDBAdaptor.fetch_slice_by_region('chromosome', '1', end=100000, strand = -1)
    #print slice.id
    #slice = coreDBAdaptor.fetch_slice_by_region('chromosome', '18', start=31129340, end=31190706)
    
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

