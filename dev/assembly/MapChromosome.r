# These functions map chromosome names from one version to another

###################################################################
# Map chromosome names
MapChromosome<-function(from, to, map) {
  # from    the original chromosome names as a vector
  # to      the names of a set of chromosomes to be mapped to
  # map     a list that maps each chromosome name to its synonyms
  
  names(map)<-toupper(names(map));
  has.it<-toupper(from) %in% names(map);
  
  out<-from; # output vector
  names(out)<-from
  out[has.it]<-sapply(from[has.it], function(fr) {
    nms<-map[toupper(fr)][[1]];
    mtch<-to[toupper(to) %in% toupper(nms)];
    if (length(mtch) == 0) fr else mtch[1];
  });
  
  dup<-out[duplicated(out)];  
  if (length(dup) > 0) out[out %in% dup]<-from[out %in% dup];
  
  out;
}

# Rename chromosome names of a GRanges object by mapping
RenameSeqlevels<-function(gr, to, map) {
  # gr      an GRanges object with chromosome names to be mapped
  # to      the names of a set of chromosomes to be mapped to
  # map     a list that maps each chromosome name to its synonyms
  
  lvl<-seqlevels(gr);
  seqlevels(gr)<-lvl[lvl %in% as.vector(seqnames(gr))];
  
  c<-MapChromosome(seqlevels(gr), to, map); 
  chr<-c[as.vector(seqnames(gr))];
  chr<-Rle(chr);
  chr@values<-as.factor(chr@values);
  gr@seqnames<-chr;
  gr@seqinfo@seqnames<-levels(chr@values);
  
  gr;
}
