# mapping a giving set of gene IDs to IDs of a different species based on NCBI homologene database
MapHomologene<-function(homolo, ids, from, to, id.type=1:4) {
  # homolo        Data frame downloaded from NCBI homologene database with 6 columns in the order of: gene cluster ID, Taxonomy ID, NCBI gene ID, Official gene symbol, GenBank ID, and RefSeq ID
  # ids           Input IDs to be mapped
  # from    Taxonomy ID of the species the input IDs belong to
  # to    Taxonomy ID of the species the input IDs to be mapped to
  # id.type   One of 4 ID types: NCBI gene, official symbol, GenBank, or RefSeq
  
  cl<-as.vector(homolo[[1]]); # gene cluster ID
  sp<-as.vector(homolo[[2]]); # species
  id<-as.vector(homolo[[id.type[1]+2]]); # gene ID
  
  from<-as.character(from);
  to<-as.character(to);
  if (!(from %in% sp)) stop('Unknown species ', from);
  if (!(to %in% sp)) stop('Unknown species ', to);
  
  x<-split(cl[id %in% ids], id[id %in% ids]);
  y<-split(id[sp %in% to], cl[sp %in% to]);
  z<-unlist(x, use.names=FALSE);
  names(z)<-rep(names(x), sapply(x, length));
  y<-y[z];
  names(y)<-names(z);
  
  out<-unlist(y, use.names=FALSE);
  names(out)<-rep(names(y), sapply(y, length));
  
  out;
}

