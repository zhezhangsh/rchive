# Use the org.Hs.eg.db package or similar package of other species to map a set of IDs to gene ID
MapOrg2Gene<-function(pkg.name, ids=c(), type) {
  # pkg.name    Name of the package
  # ids         IDs to be mapped. Use all available IDs if empty;
  # type        Type of the ID, such as 'ENSEMBL', 'REFSEQ', 'SYMBOL'. Check package manual for possible values

  if (!require(pkg.name, character.only=TRUE)) {
    source("http://bioconductor.org/biocLite.R");
    biocLite(pkg.name); 
    if(!require(pkg.name, character.only = TRUE)) stop(pkg.name, ": package not found");
  }
  
  #if (!grepl('$org', type)) fld<-paste(sub('.db$', '', pkg.name), type, sep='') else fld<-type;
  
  library(dplyr);
  fn.db<-paste(.libPaths()[1], pkg.name, 'extdata', sub('.db', '.sqlite', pkg.name), sep='/'); 
  db<-src_sqlite(fn.db);
  t<-as.data.frame(collect(tbl(db, tolower(type)))); 
  gn<-as.data.frame(collect(tbl(db, 'genes')));

  g<-gn[,2];
  names(g)<-gn[,1]; 
  mp1<-t[, 2]; 
  mp2<-as.vector(g[as.character(t[, 1])]); 
  
  #x<-eval(parse(text=fld)); 
  #mapped_genes <- mappedkeys(x);  
  #xx <- as.list(x[mapped_genes]);
  #mp1<-unlist(xx, use.names=FALSE);
  #mp2<-rep(names(xx), sapply(xx, length)); 
  mp<-split(mp2, mp1);
  
  ids<-ids[!is.na(ids)];
  ids<-ids[ids!=''];
  if (length(ids) > 0) mp<-mp[names(mp) %in% ids]; 
  
  mp<-lapply(mp, unique);
  mp<-lapply(mp, sort); 

  mp;
}