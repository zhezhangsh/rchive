# These functions map chromosome names from one version to another

# Map chromosome names
MapChromosome<-function(from, to, map) {
  # from    the original chromosome names as a vector
  # to      the names of a set of chromosomes to be mapped to
  # map     a list that maps each chromosome name to its synonyms
  
  names(map)<-toupper(names(map));
  has.it<-toupper(from) %in% names(map);
  
  out<-from; # output vector
  out[has.it]<-sapply(from[has.it], function(fr) {
    nms<-map[toupper(fr)][[1]];
    mtch<-to[toupper(to) %in% toupper(nms)];
    if (length(mtch) == 0) fr else mtch[1];
  });
  out[duplicated(out)]<-from[duplicated(out)];
  
  out;
}

# Rename chromosome names of a GRanges object by mapping
renameSeqlevels<-function(gr, to, map, len) {
  # gr      an GRanges object with chromosome names to be mapped
  # to      the names of a set of chromosomes to be mapped to
  # map     a list that maps each chromosome name to its synonyms
  
  seqlevels(gr)<-MapChromosome(seqlevels(gr), to, map); 
  
  gr;
}