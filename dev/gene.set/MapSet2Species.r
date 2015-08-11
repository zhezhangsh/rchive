# Make a list of gene sets species-specific
MapSet2Species<-function(lst, species, path.gene=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r', sep='/')) {
  # lst       A list of gene sets
  # species   Species to be mapped to
  # path.gene Folder of files providing species-specific gene annotation
  
  # Loading in genes of each species
  fn.gene<-dir(path.gene);
  fn.species<-paste(tolower(species), '_genes_full.rds', sep='');
  names(fn.species)<-tolower(species);
  fn.species<-fn.species[fn.species %in% fn.gene];
  if (length(fn.species) == 0) NA else {
    fns<-paste(path.gene, fn.species, sep='/');
  }
  ids<-lapply(fns, function(fn) rownames(readRDS(fn)));
  names(ids)<-names(fn.species);
  
  all<-unlist(lst, use.names=FALSE);
  names(all)<-rep(names(lst), sapply(lst, length));
  mp<-lapply(ids, function(id) all[all %in% id]);
  names(mp)<-names(ids);
  mp<-mp[sapply(mp, length)>0];
  if (length(mp) == 0) NA else {
    by.spe<-lapply(mp, function(mp) {
      l<-split(as.vector(mp), names(mp));
      l[names(lst[names(lst) %in% names(l)])]; 
    });
    names(by.spe)<-names(mp);

    by.spe;
  }
}