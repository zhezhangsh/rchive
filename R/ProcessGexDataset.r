# Pre-process of gene expression data set
ProcessGexDataset<-function(gex, grp2smp, id, taxid=NA, homologene=NA, map.to='9606') {
  # gex         The original normalized and logged data set
  # grp2smp     Samples grouped 
  # id          Dataset ID
  # taxid       Taxatomy ID of the original dataset; used for homolog mapping
  # homologene  Downloaded and parsed Homologene table
  # map.to      The species to map the species of original dataset
  
  library(rchive);
  
  # Add group means to matrix
  grp2smp<-lapply(grp2smp, function(smp) smp[smp %in% colnames(gex)]);
  grp2smp<-grp2smp[sapply(grp2smp, length)>0];
  if (length(grp2smp) > 0) {
    ms<-sapply(names(grp2smp), function(nm) rowMeans(gex[, grp2smp[[nm]], drop=FALSE], na.rm=TRUE));
    gex<-cbind(gex, round(ms, 6));
  }
  gex<-cbind(gex, round(rowMeans(gex, na.rm=TRUE), 6));
  colnames(gex)[ncol(gex)]<-id;
  
  pct<-apply(gex, 2, function(g) round(dplyr::percent_rank(g)*100, 4));
  
  out<-list(list(logged=gex, percentile=pct));
  names(out)<-taxid;
  
  map.to<-map.to[map.to!=taxid];
  map.to<-map.to[!is.na(map.to)];
  
  # Add homolog 
  if (!identical(homologene, NA) & !identical(taxid, NA) & length(map.to)>0) {
    homolog<-lapply(map.to, function(tid) {
      mp<-MapHomologene(homologene, rownames(gex), taxid, tid);
      g<-apply(gex[names(mp), , drop=FALSE], 2, function(g) sapply(split(g, as.vector(mp)), mean));
      p<-apply(g, 2, function(g) round(dplyr::percent_rank(g)*100, 4));
      list(logged=g, percentile=p);
    });
    names(homolog)<-map.to;
    out[[map.to]]<-homolog[[1]];
  }
  
  out;
}