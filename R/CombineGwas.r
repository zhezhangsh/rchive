# Combine more than one collections of GWAS results
CombineGwas<-function(path.out=paste(Sys.getenv('RCHIVE_HOME'), 'data/gwas/r', sep='/'), 
  paths=paste(Sys.getenv('RCHIVE_HOME'), 'data/gwas/public', c('dbgap/r', 'gwascatalog/r'), sep='/')
  ) {
  # path.out    Path to output files
  # paths       Paths to GWAS collections
     
  if (!file.exists(path.out)) dir.create(path.out, recursive=TRUE);

  library(stringr);
  
  ##################
  # depended files #
  ##################
  # required keyword files from PubTator
  fn.key<-paste(Sys.getenv("RCHIVE_HOME"), '/data/literature/public/pubtator/r/pubmed2', 
                c('gene', 'disease', 'chemical'), '.rds', sep='');
  
  # required GWAS metadata files in each collection
  fn.req<-c("analysis.rds", "study.rds", "pubmed.rds");
  fn.missed<-lapply(paths, function(path) {
    fn<-paste(path, fn.req, sep='/');
    fn[!file.exists(fn)];
  });
  if (length(unlist(fn.missed)) > 0) 
    stop("Error: missing required file(s):\n", paste(unlist(fn.missing, use.names=FALSE), collapse='\n'));
  
  #######################
  # Loading in metadata #
  #######################
  
  ################################################################################
  ################################################################################
  ### load and combine annotation tables
  # Column names of annotation tables
  cnm<-list(
    analysis=c("Source", "Name", "Study_ID", "PubMed_ID", "URL", "Num_Variants", "Min_Phred", "Max_Phred", "Description"),
    study=c("Source", "Name", "Num_Analyses", "URL", "Description"),
    pubmed=c('Title', 'Journal', 'Date', 'URL', 'Abstract')
  );
  # Load, combine, and return metadata tables
  meta<-sapply(names(cnm), function(nm) {
    tbl.set<-lapply(paths, function(path) readRDS(paste(path, '/', nm, '.rds', sep='')));
    # combined analysis tables with the same columns
    tbls<-lapply(tbl.set, function(tbl) {
      c<-setdiff(cnm[[nm]], colnames(tbl));
      if (length(c)>0) tbl<-cbind(tbl, sapply(c, function(c) rep('', nrow(tbl))));
      tbl[, cnm[[nm]]];
    });
    
    # Re-format
    tbl<-do.call('rbind', tbls);
    rnm<-unlist(lapply(tbl.set, rownames));
    tbl<-tbl[!duplicated(rnm), ];
    rownames(tbl)<-rnm[!duplicated(rnm)];
    
    for (i in 1:ncol(tbl)) tbl[[i]]<-str_trim(as.vector(tbl[[i]]));
    if ('Description' %in% cnm[[nm]]) tbl[tbl$Description=='', 'Description']<-tbl[tbl$Description=='', 'Name'];

    saveRDS(tbl, file=paste(path.out, '/', nm, '.rds', sep=''));
    
    tbl;
  }); 
  
  ################################################################################
  ################################################################################
  ### Load and combine annotation as lists
  meta.as.lst<-lapply(names(meta), function(nm) {
    # load from individual sources
    lst.set<-lapply(paths, function(path) {
      fn.lst<-paste(path, '/', nm, '_by_id.rds', sep='');
      if (file.exists(fn.lst)) readRDS(fn.lst) else { # if not exists in the source
        tbl<-readRDS(paste(path, '/', nm, '.rds', sep=''));
        by.id<-lapply(rownames(tbl), function(id) c(ID=id, tbl[id, ]));
        names(by.id)<-rownames(tbl);
        by.id
      }
    });
    # combine lists
    lst<-do.call('c', lst.set);
    saveRDS(lst, file=paste(path.out, '/', nm, '_by_id.rds', sep=''));
    lst;
  }); 
  names(meta.as.lst)<-names(meta);
  
  ################################################################################
  ################################################################################
  ## Study/analysis to PubMed
  
  map2pubmed<-lapply(c('analysis', 'study'), function(nm) {
    mp2pm<-lapply(meta.as.lst[[nm]], function(lst) 
      if('pubmed' %in% tolower(names(lst))) lst[tolower(names(lst))=='pubmed'][[1]] else '');
    
    # if available directly from source
    for (i in 1:length(paths)) {
      fn.mp<-paste(paths[[i]], '/', nm, '2pubmed.rds', sep='');
      if (file.exists(fn.mp)) {
        mp<-readRDS(fn.mp);
        mp<-mp[names(mp) %in% names(mp2pm)];
        if (length(mp) > 0) mp2pm[names(mp)]<-mp;
      };
    }
    
    saveRDS(lapply(mp2pm, function(mp) sort(unique(mp))), file=paste(path.out, '/', nm, '2pubmed.rds', sep=''));
    
    # Map pubmed ID to analysis and study
    pm<-unlist(mp2pm, use.names=FALSE);
    id<-rep(names(mp2pm), sapply(mp2pm, length));
    pm2id<-lapply(split(id, pm), function(x) sort(unique(x)));
    pm2id<-pm2id[names(pm2id) %in% rownames(readRDS(paste(path.out, 'pubmed.rds', sep='/')))];
    saveRDS(pm2id, file=paste(path.out, '/pubmed2', nm, '.rds', sep=''));    
                  
    mp2pm;
  });
   
  pm.all<-sort(unique(unlist(map2pubmed)));
  
  ############
  # Keywords #
  ############
  
  ################################################################################
  ################################################################################
  # Use PubTator to map Pubmed ID to keywords
  all.keys<-lapply(fn.key, function(fn) {
    ky<-readRDS(fn);
    ky[names(ky) %in% pm.all];
  });
  saveRDS(all.keys, file=paste(path.out, 'keyword_all.rds', sep='/'));

  ##########################################
  # Gene keywords
  gn.id<-strsplit(gn$EntrezGeneID, ',');
  n.gn.id<-sapply(gn.id, length);
  gn.id<-unlist(gn.id);
  key.gn<-data.frame(PMID=rep(gn$PMID, n.gn.id), ID=paste('ENTREZ:', gn.id, sep=''), Type='Gene', URL=paste("http://www.ncbi.nlm.nih.gov/gene/?term=", gn.id, sep=''), Synonym=rep(as.vector(gn[, 4]), n.gn.id), stringsAsFactors=FALSE);
  
  ##########################################
  # Disease keywords
  ds<-readRDS(pt.disease);
  ds<-ds[ds[,1] %in% pm.all, ];
  
  id<-sapply(strsplit(ds[,2], ':'), function(x) x[2]);
  db<-sapply(strsplit(ds[,2], ':'), function(x) x[1]);
  url<-rep('', length(id));
  url[db=='MESH']<-paste('http://www.ncbi.nlm.nih.gov/mesh/?term=', id[db=='MESH'], sep='');
  url[db=='OMIM']<-paste('http://omim.org/entry/', id[db=='OMIM'], sep='');
  key.ds<-data.frame(PMID=ds[,1], ID=ds[,2], Type='Disease', URL=url, Synonym=as.vector(ds[, 3]), stringsAsFactors=FALSE);
  
  ##########################################
  # Chemical keywords
  cm<-readRDS(pt.chemical);
  cm<-cm[cm[,1] %in% pm.all, ];
  
  id<-sapply(strsplit(cm[,2], ':'), function(x) x[2]);
  db<-sapply(strsplit(cm[,2], ':'), function(x) x[1]);
  url<-rep('', length(id));
  url[db=='MESH']<-paste('http://www.ncbi.nlm.nih.gov/mesh/?term=', id[db=='MESH'], sep='');
  url[db=='CHEBI']<-paste('http://www.ebi.ac.uk/chebi/advancedSearchFT.do?searchString=', id[db=='CHEBI'], sep='');
  key.cm<-data.frame(PMID=cm[,1], ID=cm[,2], Type='Chemical', URL=url, Synonym=as.vector(cm[, 3]), stringsAsFactors=FALSE);
  key.cm<-key.cm[key.cm[,4]!='', ];
  
  
}