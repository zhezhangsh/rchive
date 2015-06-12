# Wrap up a gene expression collection for the GeEx APP 
PrepareGeexCollection<-function(path.coll, 
                                path.gene=paste(Sys.getenv('RCHIVE_HOME'), 'data/gene/public/entrez/r', sep='/'), 
                                path.homo=paste(Sys.getenv('RCHIVE_HOME'), 'data/gene/public/homologene/r', sep='/')) {
  # path.coll   Path to the folder of the existing data collection
  
  library(rchive);
  library(awsomics);
  
  ############## Required files
  fn.meta<-paste(path.coll, 'metadata.rds', sep='/');
  fn.gex<-paste(path.coll, 'gex.rds', sep='/');
  if (!file.exists(fn.meta)) 
    stop("Error: metadata file not exist:", fn.meta, '\n');
  if (!file.exists(fn.gex)) 
    stop("Error: expression data file not exist:", fn.gex, '\n');
  ##############
  
  # Loading in data
  meta<-readRDS(fn.meta);
  gex<-readRDS(fn.gex);
  ds<-meta$metadata$Dataset;
  ds.id<-intersect(names(gex), rownames(ds));
  if (length(ds.id) == 0) stop('Error: no matching dataset ID in metadata table and expression matrix\n');
  ds<-ds[ds.id, ];
  gex<-gex[ds.id];
  tax.nm<-tolower(ds$Species);
  names(tax.nm)<-rownames(ds);
  
  ############## Required files
  fn.homo<-paste(path.homo, 'homologene.rds', sep='/');
  if (!file.exists(fn.homo)) 
    stop("Error: Homologene file not exist:", fn.homo, '\n');
  ##############
    
  ############## Required files
  fn.gene<-paste(path.gene, '/', tax.nm, '_genes_full.rds', sep='');
  fn.tax<-paste(path.gene, 'name2taxid.rds', sep='/');
  if (!file.exists(fn.tax)) 
    stop("Error: Species name to Tax ID mapping not exist:", fn.tax, '\n');
  if (length(fn.gene[!file.exists(fn.gene)]) > 0) 
    stop('Error: gene information file not exist:', paste(unique(tax.nm[!file.exists(fn.gene)]), collapse='; '));
  tax.id<-readRDS(fn.tax);
  names(tax.id)<-tolower(names(tax.id));
  x<-setdiff(tax.nm, names(tax.id));
  if (length(x) > 0) stop("Error: species not included in gene repository:", paste(x, collapse='; '), '\n');
  taxid<-tax.id[tax.nm];
  ##############
  
  ##################################################################################################
  # Dataset, group, sample metadata
  id2nm<-do.call('c', lapply(meta$metadata, function(meta) meta[, 'Name']));
  names(id2nm)<-unlist(lapply(meta$metadata, rownames));
  mapping<-list(id2name=id2nm);
  
  # Full gene annotation
  anno.all<-do.call('rbind', lapply(fn.gene, readRDS));
  gn.id<-lapply(gex, rownames);
  gn.id<-lapply(gn.id, function(id) id[id %in% rownames(anno.all)]);
  gn.id<-lapply(split(gn.id, tax.nm), function(id) unique(unlist(id)));
  anno<-anno.all[unlist(gn.id), ];
  anno<-cbind(Species=rep(names(gn.id), sapply(gn.id, length)), anno);
  anno[[1]]<-as.vector(anno[[1]]);

  # Pre-process data
  homo<-readRDS(fn.homo);
  data<-lapply(names(gex), function(nm) {
    cat("Processing dataset:", nm, '\n');
    g<-gex[[nm]];
    smp<-meta$metadata$Sample;
    snm<-setdiff(colnames(g), rownames(smp));
    if (length(snm) > 0) 
      stop("Error: samples not matching between expression data and metadata:", paste(snm, collapse='; '), '\n');
    smp<-smp[colnames(g), ];
    grp2smp<-split(rownames(smp), smp$Group);
    if (tax.nm[nm] == 'human') from<-'9606' else from<-as.vector(taxid[tax.nm[nm]]);
    d<-ProcessGexDataset(g, grp2smp, nm, from, homo, '9606');
    names(d)<-sapply(names(d), function(nm) names(taxid)[taxid==nm][1]);
    saveRDS(d, file=paste(path.coll, '/gex_', nm, '.rds', sep=''));
    d$human;
  });
  
  # Combine data sets into one matrix, using human gene ID
  rnm<-unique(unlist(lapply(data, function(d) rownames(d[[1]]))));
  cnm<-unique(unlist(lapply(data, function(d) colnames(d[[1]]))));
  logged<-pct<-matrix(NA, nr=length(rnm), nc=length(cnm), dimnames=list(rnm, cnm));
  gex.comb<-list(logged=logged, percentile=pct);
  
  for (i in 1:length(data)) {
    logged[rownames(data[[i]][[1]]), colnames(data[[i]][[1]])]<-data[[i]][[1]];
    pct[rownames(data[[i]][[2]]), colnames(data[[i]][[2]])]<-data[[i]][[2]];    
  }
  
  
  # all gene annotation
  anno<-anno[order(as.numeric(rownames(anno))), ];
  mapping$gene2name<-anno$Symbol;
  names(mapping$gene2name)<-rownames(anno);
  
  # Species mapping
  sp2id<-split(rownames(anno), anno$Species);
  mp2hs<-lapply(names(sp2id), function(sp) {
    if (sp == 'human') {
      id<-sp2id[[sp]];
      names(id)<-id;
      id;
    } else {
      MapHomologene(homo, sp2id[[sp]], taxid[sp], '9606');
    }
  });
  mp2hs<-do.call('c', mp2hs);
  mapping$species2human<-split(as.vector(mp2hs), names(mp2hs));
  mapping$human2species<-split(names(mp2hs), as.vector(mp2hs));

  homo.id<-setdiff(mp2hs, rownames(anno));
  if (length(homo.id) > 0) {
    anno0<-anno.all[homo.id, , drop=FALSE];
    anno0<-cbind(Species=rep('human', nrow(anno0)), anno0);
    anno0[[1]]<-as.vector(anno0[[1]]);
    anno<-rbind(anno, anno0);
  } 
  ##################################################################################################
  
  # Tables for content browsing 
  browse.tbl<-meta[[3]];
  names(browse.tbl)<-c('Data set', 'Group', 'Sample');
  gn.url<-awsomics::AddHref(rownames(anno), paste("http://www.ncbi.nlm.nih.gov/gene/?term=", rownames(anno), sep=''));  
  browse.tbl$Gene<-data.frame(gn.url, anno[, c(1, 2, 6, 7, 9, 8)], stringsAsFactors=FALSE);
  colnames(browse.tbl$Gene)<-c('ID', 'Species', 'Name', 'Chromosom', 'Location', 'Type', 'Title');
  
  saveRDS(gex.comb, file=paste(path.coll, 'gex_combined.rds', sep='/'));
  saveRDS(anno, file=paste(path.coll, 'gene.rds', sep='/'));
  saveRDS(browse.tbl, file=paste(path.coll, 'browse_table.rds', sep='/'));
  saveRDS(mapping, file=paste(path.coll, 'mapping.rds', sep='/'));
  
  mapping;
}