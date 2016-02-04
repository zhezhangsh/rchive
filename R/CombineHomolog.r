# Combine homolog mapping from multiple sources
CombineHomolog<-function(mapped, reorder=TRUE) {
  # reorder     reorder genes by gene ID order
  
  if (reorder) {
    mapped<-lapply(mapped, function(mp) lapply(mp, function(x) x[order(as.numeric(x))])); 
  }
  
  ids<-lapply(mapped, names); 
  ids<-sort(unique(unlist(ids, use.names=FALSE))); 
  tbl<-sapply(mapped, function(x) x[ids]); 
  rownames(tbl)<-ids;
  
  lst<-apply(tbl, 1, function(x) unlist(x, use.names=FALSE)); 
  lst<-lapply(lst, unique);
  
  
  fst<-sapply(lst, function(x) x[1]); 
  oth<-sapply(lst, function(x) paste(x[-1], collapse=';'));
  smm<-data.frame(Representative=fst, More=oth);
  
  list(summary=smm, list=lst, full=tbl); 
}