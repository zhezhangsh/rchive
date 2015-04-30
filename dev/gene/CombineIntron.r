CombineIntron<-function(tx2in, genome=NA) {
  # tx2in   A list of all transcript-to-intron collection, each element is a GRangesList of transcript ID, named by source name
  # genome  A BSgenome object of the genome version, used to get donor/acceptor sequences.

  # combine all introns
  all<-unlist(GRangesList(tx2in));
  all<-all[!duplicated(all)]; # remove duplicates (same chromosome, start, end, and strand)
  all<-all[order(as.vector(seqnames(all)), start(all))];
  names(all)<-1:length(all);

  # whether an intron is included by a source
  mtrx<-sapply(tx2in, function(t2i) countOverlaps(all, t2i, type='equal'));
  mtrx[mtrx>1]<-1;
  colnames(mtrx)<-names(tx2in);
  meta<-mtrx;
  
  if (!identical(genome, NA)) {
    # find chromosome not in genome and end position beyond chromosome length
    chr.len<-seqlengths(genome);
    chr.len<-c(chr.len, 'NA'=Inf);
    has.chr<-as.vector(seqnames(all)) %in% seqnames(genome);
    end<-end(all)+2;
    chr<-as.vector(seqnames(all));
    chr[!has.chr]<-'NA';
    len<-chr.len[chr];
    
    # has the sequence in genome
    has.seq<-has.chr & len>=end & start(all)>2
    
    # the donor/acceptor bases
    left<-GRanges(seqnames(all), IRanges(start(all)-2, start(all)-1), strand(all));
    right<-GRanges(seqnames(all), IRanges(end(all)+1, end(all)+2), strand(all));
    #seq1<-getSeq(genome, left[has.seq]);
    #seq2<-getSeq(genome, right[has.seq]);
    mtrx<-matrix('', nr=length(all), nc=4, dimnames=list(1:length(all), c('exon_left', 'intron_left', 'intron_right', 'exon_right')));
    mtrx[has.seq, 1]<-as.character(getSeq(genome, left[has.seq]));
    mtrx[has.seq, 2]<-as.character(getSeq(genome, shift(left, 2)[has.seq]));
    mtrx[has.seq, 3]<-as.character(getSeq(genome, shift(right, -2)[has.seq]));
    mtrx[has.seq, 4]<-as.character(getSeq(genome, right[has.seq]));
    
    meta<-data.frame(meta, mtrx, stringsAsFactors=FALSE);
  }
   
  elementMetadata(all)<-meta;
  
  all;
}