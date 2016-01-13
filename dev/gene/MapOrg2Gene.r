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
  
  if (!grepl('$org', type)) fld<-paste(sub('.db$', '', pkg.name), type, sep='') else fld<-type;
  
  x<-eval(parse(text=fld));
  mapped_genes <- mappedkeys(x);
  xx <- as.list(x[mapped_genes]);
  
  mp1<-unlist(xx, use.names=FALSE);
  mp2<-rep(names(xx), sapply(xx, length)); 
  mp<-split(mp2, mp1);
  
  ids<-ids[!is.na(ids)];
  ids<-ids[ids!=''];
  if (length(ids) > 0) mp<-mp[names(mp) %in% ids]; 
  
  mp<-lapply(mp, unique);
  mp<-lapply(mp, sort); 

  mp;
}