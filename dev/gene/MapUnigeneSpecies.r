# Map genes between species using parsed Unigene information
MapUnigeneSpecies<-function(mapped, pkg.name) {
  # mapped    Parsed table of Unigene mapping
  # pkg.name  org package of the species to be mapped to
  
  library(rchive);
  
  if (!require(pkg.name, character.only=TRUE)) {
    source("http://bioconductor.org/biocLite.R");
    biocLite(pkg.name); 
    if(!require(pkg.name, character.only = TRUE)) stop(pkg.name, ": package not found");
  }
  
  sp<-sub('.db', 'ORGANISM', pkg.name); 
  sp<-eval(parse(text=sp));
  sp<-as.character(genomes::ncbiTaxonomy(sp)[1,1]);
  
  mapped<-mapped[mapped[,2]==sp, , drop=FALSE]; 
  if (nrow(mapped) == 0) list() else {
    prot<-mapped$PROTID;
    nc<-nchar(prot);
    brk<-regexpr('\\.', prot); 
    brk[brk==-1]<-nc[brk==-1]+1;
    prot<-substr(prot, 1, brk-1); 
  }
  
  out<-as.list(rep('', length(prot))); 
  names(out)<-prot; 
  
  mp<-MapOrg2Gene(pkg.name, prot, 'REFSEQ'); 
  out[names(mp)]<-mp; 
  names(out)<-as.vector(mapped[, 1]); 
  
  len<-sapply(out, length); 
  ids<-rep(names(out), len); 
  gns<-unlist(out, use.names=FALSE); 
  pct<-rep(mapped$PCT, len); 
  names(pct)<-gns; 
  
  mp<-split(pct, ids); 
  mp<-lapply(mp, function(mp) sort(mp, decreasing = TRUE)); 
  mp<-mp[names(mp) != '']; 
  mp; 
}