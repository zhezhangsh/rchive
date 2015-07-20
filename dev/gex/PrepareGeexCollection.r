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
  fn.browse<-paste(path.coll, 'browse_table.rds', sep='/');
  if (!file.exists(fn.meta)) 
    stop("Error: metadata file not exist:", fn.meta, '\n');
  if (!file.exists(fn.gex)) 
    stop("Error: expression data file not exist:", fn.gex, '\n');
  if (!file.exists(fn.browse)) 
    stop("Error: formatted tables not exist:", fn.gex, '\n');
  
  ##############
  
  # Loading in data
  meta<-readRDS(fn.meta);
  gex<-readRDS(fn.gex);
  ds<-meta$Dataset;
  ds.id<-intersect(names(gex), rownames(ds));
  if (length(ds.id) == 0) stop('Error: no matching dataset ID in metadata table and expression matrix\n');
  ds<-ds[ds.id, ];
  gex<-gex[ds.id];
  gex<-lapply(gex, as.matrix);
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
  id2nm<-do.call('c', lapply(meta, function(meta) meta[, 'Name']));
  names(id2nm)<-unlist(lapply(meta, rownames));
  nm2id<-names(id2nm);
  mapping<-list(
    id2name=id2nm,
    name2id=setNames(names(id2nm), as.vector(id2nm)),
    id2longname=setNames(paste(names(id2nm), id2nm, sep=': '), names(id2nm)),
    longname2id=setNames(names(id2nm), paste(names(id2nm), id2nm, sep=': ')),
    name2longname=setNames(paste(names(id2nm), id2nm, sep=': '), as.vector(id2nm)),
    longname2name=setNames(as.vector(id2nm), paste(names(id2nm), id2nm, sep=': '))
  );
  
  # Pre-process data
  homo<-readRDS(fn.homo);
  data<-lapply(names(gex), function(nm) {
    cat("Processing dataset:", nm, '\n');
    g<-gex[[nm]];
    smp<-meta$Sample;
    snm<-setdiff(colnames(g), rownames(smp));
    if (length(snm) > 0) 
      stop("Error: samples not matching between expression data and metadata:", paste(snm, collapse='; '), '\n');
    smp<-smp[colnames(g), ];
    grp2smp<-split(rownames(smp), smp$Group);
    if (tax.nm[nm] == 'human') from<-'9606' else from<-as.vector(taxid[tax.nm[nm]]);
    d<-ProcessGexDataset(g, grp2smp, nm, from, homo, '9606');
    names(d)<-sapply(names(d), function(nm) names(taxid)[taxid==nm][1]);
    saveRDS(d, file=paste(path.coll, '/gex_', nm, '.rds', sep=''));
    lapply(1:2, function(i) do.call('rbind', lapply(d, function(d) d[[i]])));    
  });
  
  # Combine data sets into one matrix, using human gene ID
  cnm<-unique(sort(unique(unlist(lapply(data, function(d) colnames(d[[1]]))))));
  rnm<-unique(sort(as.numeric(unlist(lapply(data, function(d) rownames(d[[1]]))))));
  logged<-pct<-matrix(NA, nr=length(rnm), nc=length(cnm), dimnames=list(rnm, cnm));
  for (i in 1:length(data)) {
    logged[rownames(data[[i]][[1]]), colnames(data[[i]][[1]])]<-data[[i]][[1]];
    pct[rownames(data[[i]][[2]]), colnames(data[[i]][[2]])]<-data[[i]][[2]];    
  }
  gex.comb<-list(logged=logged, percentile=pct);
  
  # Full gene annotation
  names(fn.gene)<-tax.nm;
  fn.gene<-fn.gene[!duplicated(fn.gene)];
  anno.all<-lapply(fn.gene, readRDS);
  sp<-rep(names(fn.gene), sapply(anno.all, nrow));
  id<-as.vector(unlist(lapply(anno.all, rownames)));
  anno.all<-do.call('rbind', anno.all);
  anno<-data.frame(Species=sp, anno.all, row.names=id, stringsAsFactors=FALSE);
  gex.comb<-lapply(gex.comb, function(g) g[rownames(g) %in% rownames(anno), , drop=FALSE]);
  anno<-anno[rownames(gex.comb[[1]]), , drop=FALSE];
    
  # data set to gene mapping
  ds2gene<-lapply(names(gex), function(nm) {
    d<-readRDS(paste(path.coll, '/gex_', nm, '.rds', sep=''));
    lapply(d, function(d) rownames(d$logged));
  });
  names(ds2gene)<-names(gex);
  gene2ds<-split(rep(names(ds2gene), sapply(lapply(ds2gene, unlist), length)), unlist(ds2gene, use.names=FALSE));
  gene2ds<-lapply(gene2ds, function(ds) sort(unique(ds)));
  gene2ds<-gene2ds[order(as.numeric(names(gene2ds)))];
  mapping$dataset2gene<-ds2gene;
  mapping$gene2dataset<-gene2ds;
  mapping$gene2name<-anno$Symbol;
  names(mapping$gene2name)<-rownames(anno);
  mapping$name2gene<-setNames(names(mapping$gene2name), as.vector(mapping$gene2name));
  
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

#   homo.id<-setdiff(mp2hs, rownames(anno));
#   if (length(homo.id) > 0) {
#     anno0<-anno.all[homo.id[homo.id %in% rownames(anno.all)], , drop=FALSE];
#     anno0<-cbind(Species=rep('human', nrow(anno0)), anno0);
#     anno0[[1]]<-as.vector(anno0[[1]]);
#     anno<-rbind(anno, anno0);
#   } 
  
  # Count gene 2 dataset/sample
  anno<-anno[order(as.numeric(rownames(anno))), ];
  n.ds<-sapply(gene2ds[rownames(anno)], length);
  n.smp<-sapply(gene2ds[rownames(anno)], function(ds) sum(meta$Dataset[ds, 'Num_Sample']));
  anno<-data.frame(anno[, c('Species', 'Symbol')], Num_Dataset=n.ds, Num_Sample=n.smp, anno[, colnames(anno)[!(colnames(anno) %in% c('Species', 'Symbol'))], drop=FALSE], stringsAsFactors=FALSE);
  ##################################################################################################
  
  # Tables for content browsing 
  browse.tbl<-readRDS(fn.browse)[1:3];
  names(browse.tbl)<-c('Data set', 'Group', 'Sample');
  gn.url<-awsomics::AddHref(rownames(anno), paste("http://www.ncbi.nlm.nih.gov/gene/?term=", rownames(anno), sep=''));  
  browse.tbl$Gene<-data.frame(gn.url, anno[, c('Species', 'Symbol', 'Num_Dataset', 'Num_Sample', 'type_of_gene', 'description')], stringsAsFactors=FALSE);
  colnames(browse.tbl$Gene)<-c('ID', 'Species', 'Name', 'Num_Dataset', 'Num_Sample', 'Type', 'Title');
  
  saveRDS(gex.comb, file=paste(path.coll, 'gex_combined.rds', sep='/'));
  saveRDS(anno, file=paste(path.coll, 'gene.rds', sep='/'));
  saveRDS(browse.tbl, file=paste(path.coll, 'browse_table.rds', sep='/'));
  saveRDS(mapping, file=paste(path.coll, 'mapping.rds', sep='/'));
  
  mapping;
}