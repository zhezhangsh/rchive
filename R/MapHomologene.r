# Convert a data matrix of one species to another species through homologene
MapHomologeneMatrix<-function(mtrx, ids=rownames(mtrx), homolo, from, to, id.type=1:4, 
                              dup.method=c('mean', 'max', 'min', 'random')) {
  # mtrx          The original data matrix
  # homolo        Data frame downloaded from NCBI homologene database with 6 columns in the order of: gene cluster ID, Taxonomy ID, NCBI gene ID, Official gene symbol, GenBank ID, and RefSeq ID
  # ids           Input IDs to be mapped; the original gene ID. Length must match number of rows in matrix
  # from          Taxonomy ID of the species the input IDs belong to
  # to            Taxonomy ID of the species the input IDs to be mapped to
  # id.type       One of 4 ID types: NCBI gene, official symbol, GenBank, or RefSeq
  # dup.method    What to do with multiple rows mapped to the same gene
  
  if (length(ids) != nrow(mtrx)) stop('Error: length of IDs not matching number of rows\n');
  mtrx<-as.matrix(mtrx);
  rownames(mtrx)<-1:nrow(mtrx);
  
  # Map IDs
  mp<-rchive::MapHomologene(homolo, ids, from, to, id.type);
  mp<-mp[ids];
  names(mp)<-rownames(mtrx);
  id2row<-split(names(mp), mp);
  
  # Generate the new matrix
  dup.method<-tolower(dup.method)[1];
  d<-sapply(id2row, function(rid) {
    if (length(rid) == 1) mtrx[rid, ] else {
      if (dup.method=='random') mtrx[sample(rid, 1), ] else 
        if (dup.method=='min') mtrx[rid[order(rowMeans(mtrx[rid, ], na.rm=TRUE))][1], ] else 
          if (dup.method=='min') mtrx[rid[order(rowMeans(mtrx[rid, ], na.rm=TRUE))][length(rid)], ] else 
            colMeans(mtrx[rid, ], na.rm=TRUE);
    }
  });
  d<-t(d);
  
  h<-homolo[homolo[, 'Taxonomy_ID']==to, 3:6];
  rownames(h)<-h[, id.type];
  anno<-h[rownames(d), , drop=FALSE];
  
  list(anno=anno, data=d);
}

# mapping a giving set of gene IDs to IDs of a different species based on NCBI homologene database
MapHomologene<-function(homolo, ids=c(), from, to, id.type=1:4) {
  # homolo        Data frame downloaded from NCBI homologene database with 6 columns in the order of: gene cluster ID, Taxonomy ID, NCBI gene ID, Official gene symbol, GenBank ID, and RefSeq ID
  # ids           Input IDs to be mapped
  # from          Taxonomy ID of the species the input IDs belong to
  # to            Taxonomy ID of the species the input IDs to be mapped to
  # id.type       One of 4 ID types: NCBI gene, official symbol, GenBank, or RefSeq
  
  cl<-as.vector(homolo[[1]]); # gene cluster ID
  sp<-as.vector(homolo[[2]]); # species
  id<-as.vector(homolo[[id.type[1]+2]]); # gene ID
  
  from<-as.character(from);
  to<-as.character(to);
  if (from=='559292') from<-'4932';
  if (to=='559292') to<-'4932';  
  if (!(from %in% sp)) stop('Unknown species ', from);
  if (!(to %in% sp)) stop('Unknown species ', to);
  
  if (length(ids) > 0) x<-split(cl[id %in% ids], id[id %in% ids]) else x<-split(cl, id);
  y<-split(id[sp %in% to], cl[sp %in% to]);
  z<-unlist(x, use.names=FALSE);
  names(z)<-rep(names(x), sapply(x, length));
  y<-y[z];
  names(y)<-names(z);
  
  out<-unlist(y, use.names=FALSE);
  names(out)<-rep(names(y), sapply(y, length));
  
  out;
}

