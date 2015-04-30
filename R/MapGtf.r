# Mapping between types of regions based on parsed GTF file
GtfTranscript2Exon<-function(gr, cnm=c('parent', 'transcript_id'), cnm.id=c('id', 'transcript_id'), intron=TRUE) {
  # gr      A GTF file parsed into a GRanges object
  # cnm     Column name of the ID of transcript an exon is mapped to. NCBI uses <parent> for column name.
  # cnm.id  Column name of the unique ID of a transcript. NCBI uses <id> for column name
  #         It's known that NCBI sometimes mistakenly gives the same tRNAs multiple lines so the <id> column is not unique
  # intron  If TRUE, get introns based on exon locations. Note that it needs not only exon locations but also transcript-to-exon mapping to get introns.
  cnm<-cnm[cnm %in% colnames(elementMetadata(gr))];
  cnm.id<-cnm.id[cnm.id %in% colnames(elementMetadata(gr))];
  
  ex<-gr[gr$type %in% 'exon']; # all exons
  
  if (length(ex)==0 | length(cnm)==0 | length(cnm.id)==0) {
    warning('Metadata columns not found');
    NA;
  } else {
    ex.id<-names(ex); # row ID of exons
    pr<-as.vector(elementMetadata(ex)[[cnm[1]]]);
    
    # if more than one parent per exon (unlikely to happen)
    if (is.list(pr)) {
      ex.id<-rep(ex.id, sapply(pr, length));
      ex<-ex[ex.id];
      pr<-unlist(pr, use.names=FALSE);
    }
    
    # Check ambiguity
    chr<-as.vector(seqnames(ex));
    str<-as.vector(strand(ex));
    pr2chr<-lapply(split(chr, pr), function(x) unique(x));
    pr2str<-lapply(split(str, pr), function(x) unique(x));
    n1<-sapply(pr2chr, length);
    n2<-sapply(pr2str, length);
    if (max(n1)>1) warning('Ambiguous mapping: exons on different chromosomes have the same transcript ID.');
    if (max(n2)>1) warning('Ambiguous mapping: exons on different strands have the same transcript ID.')
    
    # annotation info. of transcripts themselves
    tx<-gr[as.vector(elementMetadata(gr)[[cnm.id[1]]]) %in% unique(pr) & gr$type!='exon'];
    tx.id<-as.vector(elementMetadata(tx)[[cnm.id[1]]]);
    
    if (length(setdiff(pr, tx.id))<0.1*length(pr) & length(tx.id)<0.9*length(pr)) { # most transcripts have their own entries in the GTF file
      # mapping 
      if (max(n1)==1 & max(n2)==1) { # if no ambiguity
        pr.id<-names(tx);
        names(pr.id)<-as.vector(elementMetadata(tx)[[cnm.id[1]]]);
        pr.id<-pr.id[pr];
      } else {
        olap<-as.matrix(findOverlaps(ex, tx, type='within'));
        olap<-olap[elementMetadata(ex)[[cnm[1]]][olap[,1]] == elementMetadata(tx)[[cnm.id[1]]][olap[,2]], ];
        pr.id<-names(tx)[olap[,2]];
        names(pr.id)<-as.vector(elementMetadata(tx)[[cnm.id[[1]]]])[olap[,2]];
        ex<-ex[olap[,1]];    
      }
      
      # preparing outputs: transcript2exon mapping and annotation of transcripts
      elementMetadata(ex)<-NULL;
      pr2ex<-split(ex, pr.id);
      pr2stt<-split(start(ex), pr.id); # Specify transcript location based on its exons
      pr2end<-split(end(ex), pr.id);
      pr2wid<-split(width(ex), pr.id);
      pr2chr<-sapply(split(as.vector(seqnames(ex)), pr.id), unique);
      pr2str<-sapply(split(as.vector(strand(ex)), pr.id), unique);
      stt<-sapply(pr2stt, min);
      end<-sapply(pr2end, max);
      wid<-sapply(pr2wid, sum);
      anno<-GRanges(unlist(pr2chr), IRanges(stt, end), strand=unlist(pr2str), length=wid, n_exon=elementLengths(pr2ex));
      anno<-anno[order(as.vector(seqnames(anno)), start(anno))];
      elementMetadata(anno)<-cbind(elementMetadata(anno), elementMetadata(tx[names(anno)]));      
    } else { # just entries for transcripts
      pr.id<-paste(as.vector(elementMetadata(ex)[[cnm.id[1]]]), as.vector(seqnames(ex)), as.integer(strand(ex))-1, sep='_')
      tx.id<-as.vector(elementMetadata(ex)[[cnm.id[1]]]);
      elementMetadata(ex)<-NULL;
      pr2ex<-split(ex, pr.id);
      pr2stt<-split(start(ex), pr.id); # Specify transcript location based on its exons
      pr2end<-split(end(ex), pr.id);
      pr2wid<-split(width(ex), pr.id);
      pr2chr<-sapply(split(as.vector(seqnames(ex)), pr.id), unique);
      pr2str<-sapply(split(as.vector(strand(ex)), pr.id), unique);
      pr2tx.id<-sapply(split(tx.id, pr.id), function(x) x[1]);
      stt<-sapply(pr2stt, min);
      end<-sapply(pr2end, max);
      wid<-sapply(pr2wid, sum);
      anno<-GRanges(unlist(pr2chr), IRanges(stt, end), strand=unlist(pr2str),  transcript_id=pr2tx.id[names(wid)], length=wid, n_exon=elementLengths(pr2ex));
      anno<-anno[order(as.vector(seqnames(anno)), start(anno))];
      pr2ex<-pr2ex[names(anno)];
    }
    
    # output
    out<-list(transcript=anno, transcript2exon=pr2ex[names(anno)]);
    
    if (intron) {
      out$transcript2intron=GtfExon2Intron(pr2ex[names(anno)]);
    }
    
    out;
  } # end of else
}

# Get intron locations based on exon locations and transcript-to-exon mapping
GtfExon2Intron<-function(tx2ex) {
  # tx2ex   A list of transcript to exon mapping named with transcript ID. Each element is a GRanges of exons
  
  # remove single-exon transcripts (no intron);
  tx2ex<-tx2ex[elementLengths(tx2ex)>1];
  
  if (length(tx2ex) == 0) {
    NA;
  } else {
    nm0<-names(tx2ex);
    names(tx2ex)<-1:length(tx2ex);
    names(nm0)<-names(tx2ex);
    
    # merge all exons to speed up
    ex<-unlist(tx2ex, use.names=FALSE);
    ex$tx<-rep(names(tx2ex), elementLengths(tx2ex));
    ex<-ex[order(ex$tx, start(ex))]; # sort by transcript first, then location
    
    stt<-end(ex)[-length(ex)]+1; # intron start position
    end<-start(ex)[-1]-1; # intron end position
    same.tx<-ex$tx[-1] == ex$tx[-length(ex)]; # whether 2 exons belong to the same transcript
    
    # remove rows not belong to the same transcript
    intron<-ex[-1][same.tx];
    stt<-stt[same.tx];
    end<-end[same.tx];
    
    # remove conflict exons
    wid<-end-stt+1;
    intron<-intron[wid>0];
    stt<-stt[wid>0];
    end<-end[wid>0];
    
    end(intron)<-end;
    start(intron)<-stt;
    
    # transcript to intron mapping
    tx.id<-as.vector(intron$tx);
    elementMetadata(intron)<-NULL;
    tx2in<-split(intron, tx.id);
    names(tx2in)<-as.vector(nm0[names(tx2in)]); 
    
    tx2in;
  }
  
}