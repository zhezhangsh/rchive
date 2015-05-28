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
  fn.key<-paste(Sys.getenv("RCHIVE_HOME"), '/data/literature/public/pubtator/r/original_pubmed2', 
                c('gene', 'disease', 'chemical'), '.rds', sep='');
  names(fn.key)<-c('Gene', 'Disease', 'Chemical');
  
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
    ky[ky[,1] %in% pm.all, ];
  });
  names(all.keys)<-names(fn.key);
  saveRDS(all.keys, file=paste(path.out, 'keyword_all.rds', sep='/'));

  # combine keywords into a single table
  tbl<-cbind(do.call('rbind', all.keys), Type=rep(names(fn.key), sapply(all.keys, nrow)));
  tbl<-cbind(tbl[, 1:2], Type=rep(names(fn.key), sapply(all.keys, nrow)), tbl[, 3:ncol(tbl)]);

  # expand keyword IDs if multiple IDs in the same cell
  key.ids<-strsplit(tbl$MeshID, '[/,]'); # All keyword IDs
  key.rnm<-rep(rownames(tbl), sapply(key.ids, length)); # table row names
  key.ids<-unlist(key.ids, use.names=FALSE);
  key.typ<-tbl[key.rnm, 'Type']; # Keyword type
  key.pub<-tbl[key.rnm, 'PMID']; # Corresponding Pubmed ID
  key.syn<-tbl[key.rnm, 'Mentions'];
  
  # Fix keyword format as source:id
  key.ids<-sub('GeneID', 'ENTREZ', key.ids);
  key.ids[key.typ=='Gene' & !grepl(':', key.ids)]<-paste('ENTREZ:', key.ids[key.typ=='Gene' & !grepl(':', key.ids)], sep='')
  key.ids[key.typ=='Chemical' & !grepl(':', key.ids)]<-paste('MESH:', key.ids[key.typ=='Chemical' & !grepl(':', key.ids)], sep='')
  key.src<-sapply(strsplit(key.ids, ':'), function(x) x[1]);
  key.val<-sapply(strsplit(key.ids, ':'), function(x) x[2]);
  
  # Create URL
  url.pre<-c(
    'ENTREZ' = "http://www.ncbi.nlm.nih.gov/gene/?term=",
    'MESH' = "http://www.ncbi.nlm.nih.gov/mesh/?term=", 
    'OMIM' = "http://omim.org/entry/",
    'CHEBI' = "http://www.ebi.ac.uk/chebi/advancedSearchFT.do?searchString="
  )[key.src];
  url<-paste(url.pre, key.val, sep='');
  url[is.na(url.pre)]<-'';
  
  keyword<-data.frame(PMID=key.pub, ID=key.ids, Type=key.typ, URL=url, Synonym=key.syn, stringsAsFactors=FALSE);
  keys<-unique(keyword$ID);
  
  ###################################################################
  # Keyword names
  ttl<-rep('', length(keys));
  names(ttl)<-keys;
  id<-sapply(strsplit(keys,':'), function(x) x[2]);
  names(id)<-keys;
  url<-as.vector(keyword$URL);
  names(url)<-as.vector(keyword$ID);
  url<-url[keys];
  
  # Names of known genes
  f<-dir(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r/', sep='/'));
  f<-f[grep('_genes_full.rds$', f)];
  f<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r', f, sep='/')
  anno<-lapply(f, readRDS);
  anno<-do.call('rbind', anno);
  sym<-as.vector(anno[id[grep('^ENTREZ', keys)], 'Symbol']);
  desc<-as.vector(anno[id[grep('^ENTREZ', keys)], 'description']);
  nm<-paste(sym, desc, sep='; ');
  names(nm)<-keys[grep('^ENTREZ', keys)];
  nm<-nm[!is.na(sym) & !is.na(desc)];
  ttl[names(nm)]<-as.vector(nm);
  
  # Names of other keywords, parse from html pages
  library(stringr); # str_trim function
  url0<-url[names(ttl[ttl=='' & url!=''])];
  if (length(url0) > 0) {
    pgs<-GetSetURL(url0);
    names(pgs)<-names(url0);
    ttl0<-sapply(names(url0), function(nm) {
      print(nm);
      u<-url0[nm];
      #html<-strsplit(getURL(u), '\n')[[1]];
      download.file(u, 'temp.html');
      html<-scan('temp.html', what='', sep='\n', flush=TRUE);
      if (grepl('^ENTREZ:', nm)) {
        t<-strsplit(html[grep('<title>', html)], '[<[>]')[[1]][3];
        t<-sub('No items found - Gene - NCBI', '', t);
      } else if (grepl('^OMIM:', nm)) {
        mk<-paste(' ', sub('^OMIM:', '', nm), ' - ', sep='');
        t<-html[grep(mk, html)][1];
        t<-sapply(strsplit(t, mk), function(x) x[2]);
      } else if (grepl('^MESH:', nm)) {
        t<-strsplit(html[grep('<title>', html)], '[<[>]')[[1]][3];
        t<-sub(' - MeSH - NCBI', '', t);
      } else if (grepl('^CHEBI:', nm)) {
        mk<-paste('\\(', nm, '\\)$', sep='');
        t<-html[grep(mk, html)][1];
        t<-sub(mk, '', t);
      }
      str_trim(t); # take out spaces at the beginning or end
    });
    
    # If name not available, use source names
    ttl[names(ttl0)]<-ttl0;
    ttl[ttl=='No items found']<-'';
    ttl1<-ttl[is.na(ttl) | ttl==''];
    ttl.k<-keyword[keyword$ID %in% names(ttl1), ];
    mp<-split(as.vector(ttl.k$Synonym), as.vector(ttl.k$ID))[names(ttl1)];
    mp<-sapply(mp, function(mp) unique(mp)[1]);
    ttl[names(ttl1)]<-mp[names(ttl1)]; 
  }
}