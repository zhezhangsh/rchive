# Prepare data for the Meta-GWAS APP
# Combine more than one collections of GWAS results

#########################
# Depend on a few previously parsed files (see below for details)
#########################
PrepareMetaGwas<-function(phred_tbls,
                      path.out=paste(Sys.getenv('RCHIVE_HOME'), 'data/gwas/r', sep='/'), 
                      paths=paste(Sys.getenv('RCHIVE_HOME'), 'data/gwas/public', c('dbgap/r', 'gwascatalog/r'), sep='/')
  ) {
  
  # path.out    Path to output files
  # paths       Paths to GWAS collections
     
  if (!file.exists(path.out)) dir.create(path.out, recursive=TRUE);

  library(stringr);
  
  ##################
  # required files #
  ##################
  # PubMed to keyword mapping based on PubTator
  fn.key<-paste(Sys.getenv("RCHIVE_HOME"), '/data/literature/public/pubtator/r/original_pubmed2', 
                c('gene', 'disease', 'chemical'), '.rds', sep='');
  names(fn.key)<-c('Gene', 'Disease', 'Chemical');
  if (length(fn.key[!file.exists(fn.key)]) > 0) 
    stop("Error: missing depended PubTator files:\n", paste(fn.key[!file.exists(fn.key)], collapse='; ')); 
  
  ### Keyword annotation
  # Gene annotation
  fn.anno<-c(
    gene=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r/human_genes_full.rds', sep='/'),
    # Disease
    mesh.d=paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/mesh/r/mesh.rds', sep='/'),
    omim=paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/omim/r/omim.rds', sep='/'),
    # chemical
    mesh.c=paste(Sys.getenv("RCHIVE_HOME"), 'data/chemical/public/mesh/r/mesh.rds', sep='/'),
    chebi=paste(Sys.getenv("RCHIVE_HOME"), 'data/chemical/public/chebi/r/chebi.rds', sep='/')
  );
  
  if (length(fn.anno[!file.exists(fn.anno)]) > 0) 
    stop("Error: missing depended annotation files:\n", paste(fn.anno[!file.exists(fn.anno)], collapse='; ')); 
  fn.gn<-dir(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r/', sep='/'));
  fn.gn<-fn.gn[grep('_genes_full.rds$', fn.gn)];
  fn.gn<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r', fn.gn, sep='/')
  fn.chebi.alt<-paste(Sys.getenv("RCHIVE_HOME"), 'data/chemical/public/chebi/r/chebi_alternative_id.rds', sep='/')
  
  # required GWAS metadata files in each collection
  fn.req<-c("analysis.rds", "study.rds", "pubmed.rds");
  fn.missed<-lapply(paths, function(path) {
    fn<-paste(path, fn.req, sep='/');
    fn[!file.exists(fn)];
  });
  if (length(unlist(fn.missed)) > 0) 
    stop("Error: missing required file(s):\n", paste(unlist(fn.missing, use.names=FALSE), collapse='\n'));
  
  # Tables of phred scores
  fn.phred<-lapply(paths, function(path) {
    fn<-dir(path, recursive=TRUE);
    fn<-fn[grep('^phred_', fn)];
    fn<-fn[grep('.rds$', fn)];
    f<-paste(path, fn, sep='/');
    names(f)<-sub('.rds$', '', sub('^phred_', '', fn));
    f;
  });
  fn.phred<-do.call('c', fn.phred);
  fn.phred[names(phred_tbls)];
  if (length(fn.phred[is.na(fn.phred)]) > 0) 
    warning("Warning: missing phred score table(s)", paste(fn.phred[is.na(fn.phred)], collapse='\n'));
  fn.phred<-fn.phred[!is.na(fn.phred)];
  if (length(fn.phred) == 0) stop("Error: no table of phred score exists!\n");
  ################################################################################
  
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
    anno<-readRDS(paste(path.out, '/', nm, '.rds', sep=''));
    if ('PubMed_ID' %in% colnames(anno)) {
      pmid<-sapply(mp2pm, function(x) x[1]);
      names(pmid)<-names(mp2pm);
      pmid<-pmid[rownames(anno)];
      pmid[is.na(pmid)]<-'';
      anno[, 'PubMed_ID']<-pmid;
    }
    saveRDS(anno, file=paste(path.out, '/', nm, '.rds', sep=''));
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
  # Get keyword names
  anno.gn<-lapply(unique(c(fn.anno['gene'], fn.gn)), readRDS);
  id2nm.gn<-unlist(lapply(anno, function(x) x[[1]]), use.names=FALSE);
  names(id2nm.gn)<-paste('ENTREZ', unlist(lapply(anno.gn, rownames), use.names=FALSE), sep=':');
  #  
  chebi<-readRDS(fn.anno['chebi']);
  id2nm.ch<-as.vector(chebi[[1]]);
  names(id2nm.ch)<-rownames(chebi);
  chebi.alt<-readRDS(fn.chebi.alt);
  alt<-chebi[chebi.alt, 1];
  names(alt)<-names(chebi.alt);
  id2nm.ch<-c(id2nm.ch, alt[!is.na(alt)]);
  # 
  omim<-readRDS(fn.anno['omim']);
  id2nm.om<-as.vector(omim[[1]]);
  names(id2nm.om)<-paste("OMIM", rownames(omim), sep=':');
  #
  mesh<-lapply(fn.anno[c('mesh.d', 'mesh.c')], readRDS);
  id2nm.ms<-unlist(lapply(mesh, function(x) as.vector(x[[1]])), use.names=FALSE);
  names(id2nm.ms)<-paste('MESH', unlist(lapply(mesh, rownames), use.names=FALSE), sep=':');
  
  # keyword names
  id2nm<-c(id2nm.ms, id2nm.om, id2nm.ch, id2nm.gn);
  key2nm<-id2nm[keys];
  names(key2nm)<-keys;
  syn<-keyword$Synonym;
  names(syn)<-keyword$ID;
  key2nm[is.na(key2nm)]<-syn[names(key2nm[is.na(key2nm)])];
  
  ##################################################################################
  # summarize keywords
  pm2id<-split(keyword$ID, keyword$PMID);
  
  # Create key to analysis matrix
  ana<-readRDS(paste(path.out,  'analysis.rds', sep='/'));
  ana.by.id<-readRDS(paste(path.out, 'analysis_by_id.rds', sep='/'));
  ana2pm<-lapply(ana.by.id, function(x) x$PubMed);
  ana2key<-lapply(ana2pm, function(pm) unlist(pm2id[pm], use.names=FALSE))
  ana2key<-lapply(ana2key, unique)[rownames(ana)];
  key2ana<-matrix(data=0, nr=length(keys), nc=length(ana2key), dimnames=list(keys, names(ana2key)))
  for (i in 1:length(ana2key)) key2ana[ana2key[[i]], i]<-1;
  saveRDS(key2ana, file=paste(path.out, 'keyword2analysis.rds', sep='/'));

  # Create key to study matrix
  std<-readRDS(paste(path.out, 'study.rds', sep='/'));
  std.by.id<-readRDS(path.out, 'study_by_id.rds', sep='/');
  std2pm<-lapply(std.by.id, function(x) x$PubMed);
  std2key<-lapply(std2pm, function(pm) unlist(pm2id[pm], use.names=FALSE))
  std2key<-lapply(std2key, unique)[rownames(std)];
  key2std<-matrix(data=0, nr=length(keys), nc=length(std2key), dimnames=list(keys, names(std2key)))
  for (i in 1:length(std2key)) key2std[std2key[[i]], i]<-1;
  saveRDS(key2std, file=paste(path.out, 'keyword2study.rds', sep='/'));
  
  # map key to synonyms
  key2syn<-split(as.vector(keyword$Synonym), as.vector(keyword$ID))[keys];
  key2syn<-lapply(key2syn, function(x) unlist(strsplit(x, '\\|'), use.names=FALSE));
  key2syn<-lapply(key2syn, function(x) sub('^ ', '', sub(' $', '', x)));
  key2syn<-lapply(key2syn, function(x) sort(x[!duplicated(tolower(x))]));
  saveRDS(key2syn, file=paste(path.out, 'keyword2synonym.rds', sep='/'));
  
  # map key to pubmed
  key2pm<-split(keyword$PMID, keyword$ID)[keys];
  key2pm<-lapply(key2pm, unique);
  key2pm<-lapply(key2pm, sort);
  key2pub<-matrix(data=0, nr=length(key2pm), nc=length(pm.all), dimnames=list(names(key2pm), pm.all));
  for (i in 1:length(key2pm)) key2pub[names(key2pm)[i], key2pm[[i]]]<-1;
  saveRDS(key2pub, file=paste(path.out, 'keyword2pubmed.rds', sep='/'));
  
  pub2key<-apply(key2pub, 2, function(x) if (sum(x)>0)rownames(key2pub)[x>0] else c())[pm.all];
  saveRDS(pub2key, file=paste(path.out, 'pubmed2keyword.rds', sep='/'));
  
  # Keyword metadata table
  typ<-sapply(split(keyword$Type, keyword$ID), unique)[keys];
  url<-sapply(split(keyword$URL, keyword$ID), unique)[keys];
  kywd<-data.frame(Keyword=keys, Type=typ, Num_Analyses=rowSums(key2ana), Num_Study=rowSums(key2std), Num_PubMed=sapply(key2pm, length),
                   URL=url, Name=key2nm[keys], stringsAsFactors=FALSE);
  saveRDS(kywd[, -1], file=paste(path.out, 'keyword.rds', sep='/'));
  
  ######################################################################
  #### Counting
  # Count PubMed mapping
  mp2pm<-list(
    ana2pm=readRDS(file=paste(path.out, 'pubmed2analysis.rds', sep='/')),
    std2pm=readRDS(file=paste(path.out, 'pubmed2study.rds', sep='/')),
    key2pm=readRDS(file=paste(path.out, 'pubmed2keyword.rds', sep='/'))
  );
  pubmed<-readRDS(file=paste(path.out, 'pubmed.rds', sep='/'));
  ct<-sapply(mp2pm, function(x) sapply(x[rownames(pubmed)], length));
  colnames(ct)<-c('Num_Analysis', 'Num_Study', 'Num_Keyword');
  pubmed<-cbind(pubmed, ct);
  saveRDS(pubmed, paste(path.out, 'pubmed.rds', sep='/'));
  
  # Summary tables of phred collections
  ana.id<-lapply(fn.phred, function(fn) unique(colnames(readRDS(fn))));
  ana<-readRDS(paste(path.out,  'analysis.rds', sep='/'));
  std.id<-lapply(ana.id, function(x) unique(as.vector(ana[x, "Study_ID"]))); 
  std2pm<-readRDS(file=paste(path.out, 'study2pubmed.rds', sep='/'));
  pub.id<-lapply(std.id, function(x) unique(unlist(std2pm[x], use.names=FALSE)));
  tbl<-sapply(list(ana.id, std.id, pub.id), function(x) sapply(x, length));
  tbl<-data.frame(tbl, phred_tbls[names(fn.phred)], stringsAsFactors=TRUE);
  names(tbl)<-c("Num_Analysis", "Num_Study", "Num_PubMed", "Description");
  saveRDS(tbl, file=paste(path.out, 'table.rds', sep='/'));
  
  
  ################################################################################
  ################################################################################
  # create formatted DataTable for Shiny rendering
  url.tag<-function(v, url) paste('<a href="', url, '" target="_blank">', v,'</a>', sep='');
  
  tbl<-readRDS(file=paste(path.out, 'table.rds', sep='/'));
  ana<-readRDS(file=paste(path.out, 'analysis.rds', sep='/'));
  std<-readRDS(file=paste(path.out, 'study.rds', sep='/'));
  pub<-readRDS(file=paste(path.out, 'pubmed.rds', sep='/'));
  key<-readRDS(file=paste(path.out, 'keyword.rds', sep='/'));
  
  browse.tbl<-list(
    'Database tables' = cbind(ID = rownames(tbl), tbl),
    'Studies' = cbind(ID = url.tag(rownames(std), as.vector(std$URL)), std[, c('Num_Analyses', 'Name')]),
    'Analyses' = cbind(ID = url.tag(rownames(ana), as.vector(ana$URL)), ana[, c('Study_ID', 'Num_Variants', 'Min_Phred', 'Max_Phred', 'Name')]),
    'PubMed' = cbind(ID = url.tag(rownames(pub), as.vector(pub$URL)), pub[, c('Num_Analysis', 'Num_Study', 'Num_Keyword', 'Title')]),
    'Keyword' = cbind(ID = url.tag(rownames(key), as.vector(key$URL)), key[, c("Type", "Num_Analyses", "Num_Study", "Num_PubMed", "Name")])
  );
  
  saveRDS(browse.tbl, file=paste(path.out, 'browse_table.rds', sep='/'));
  
  list(Analysis=ana, Study=std, PubMed=pub, Keyword=kywd);
}