# Parse OMIM database
ParseOmim<-function(
  url=c(
    'mimTitles' = 'https://data.omim.org/downloads/2XUnXL5aSL6omABJrtJVtg/mimTitles.txt',
    'mim2gene' = 'https://omim.org/static/omim/data/mim2gene.txt',
    'morbidmap' = 'https://data.omim.org/downloads/2XUnXL5aSL6omABJrtJVtg/morbidmap.txt',
    'genemap' = 'https://data.omim.org/downloads/2XUnXL5aSL6omABJrtJVtg/genemap.txt',
    'genemap2' = 'https://data.omim.org/downloads/2XUnXL5aSL6omABJrtJVtg/genemap2.txt',
    'phenotypeseries.html' = 'https://www.omim.org/phenotypicSeriesTitle/all'
  ), 
  path.gene=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r/human_genes_synonyms2id.rds', sep='/'),
  path=paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/omim', sep='/'), 
  download.all=FALSE) {

  require(RoCA); 
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
 
  # Download from OMIM ftp site
  fn.loc <- names(url);
  fn.loc <- paste(path, 'src', fn.loc, sep='/');
  for (i in 1:length(url)) if (download.all | !file.exists(fn.loc[i])) download.file(url[i], fn.loc[i]);
  
  #######################################################################################
  # mimTitles
  mt <- scan(fn.loc[1], sep='\n', flush=TRUE, what='');
  mt <- mt[substr(mt, 1, 1)!='#'];
  mt <- t(sapply(strsplit(mt, '\t'), function(x) x[1:6]));
  rownames(mt) <- mt[, 2];
  colnames(mt) <- c('Prefix', 'Mim Number', 'Preferred_Title', 'Alternative_Title', 'Symbol_Title', 'Symbol');
  mt <- data.frame(mt[, -2], stringsAsFactors = FALSE); 

  #######################################################################################
  # mim2gene
  mm <- scan(fn.loc[2], sep='\n', flush=TRUE, what='');
  hd <- mm[grep('^#', mm)];
  hd <- hd[length(hd)];
  hd <- strsplit(sub('^# ', '', hd), '\t')[[1]];
  hd <- sapply(strsplit(hd, '\\('), function(x) x[1]);
  mm <- mm[substr(mm, 1, 1)!='#'];
  mm <- t(sapply(strsplit(mm, '\t'), function(x) x[1:5]));
  colnames(mm) <- hd; 
  rownames(mm) <- mm[, 1];
  mm <- data.frame(mm[, -1], stringsAsFactors = FALSE); 
  an <- cbind(mm, mt); 
  an$URL <- paste('https://www.omim.org/entry', rownames(an), sep='/');
  colnames(an) <- gsub('\\.', '_', colnames(an));
  colnames(an) <- sub('_$', '', colnames(an));
  saveRDS(an, file=paste(path, 'r/omim_anno.rds', sep='/'));
  
  #######################################################################################
  # morbidmap
  mb <- scan(fn.loc[3], sep='\n', flush=TRUE, what='');
  hd <- mb[grep('^#', mb)];
  hd <- hd[4];
  hd <- strsplit(sub('^# ', '', hd), '\t')[[1]];
  mb <- mb[substr(mb, 1, 1)!='#'];
  mb <- t(sapply(strsplit(mb, '\t'), function(x) x[1:length(hd)]));
  colnames(mb) <- hd; 
  mb <- data.frame(mb, stringsAsFactors = FALSE); 
  gn <- readRDS(path.gene);
  mp.id <- mb[, 3];
  mp.gn <- lapply(strsplit(mb[, 2], ', '), function(x) as.vector(gn[x[[1]]]));
  id2gn <- lapply(split(mp.gn, mp.id), function(x) unique(unlist(x))); 
  id2gn <- lapply(id2gn, function(x) x[x!='']);
  id2gn0 <- id2gn[sapply(id2gn, length)>0]; 
  gn2id <- split(rep(names(id2gn), sapply(id2gn, length)), unlist(id2gn));
  gn2id0 <- lapply(gn2id, unique); 
  saveRDS(mb, file=paste(path, 'r/morbid_anno.rds', sep='/'));
  saveRDS(id2gn0, file=paste(path, 'r/morbid_omim2gene.rds', sep='/'));
  saveRDS(gn2id0, file=paste(path, 'r/morbid_gene2omim.rds', sep='/'));
  
  #######################################################################################
  # genemap
  mp <- scan(fn.loc[4], sep='\n', flush=TRUE, what='');
  hd <- strsplit(sub('^#', '', mp[4]), '\t')[[1]]; 
  mp <- mp[substr(mp, 1, 1)!='#'];
  mp <- strsplit(mp, '\t');
  mp <- sapply(mp, function(x) x[1:length(hd)]); 
  mp <- t(mp); 
  colnames(mp) <- hd; 
  mp.gn<-strsplit(mp[, 'Gene Symbols'], ', ');
  mp.id<-mp[, 'MIM Number'];
  mp.id<-rep(mp.id, sapply(mp.gn, length));
  mp.gn<-unlist(mp.gn, use.names=FALSE);
  hu.gn<-readRDS(path.gene);
  gn<-hu.gn[mp.gn];
  mp.id<-rep(mp.id, sapply(gn, length));
  mp.gn<-unlist(gn, use.names=FALSE);
  id2gn <- lapply(split(mp.gn[mp.gn!=''&mp.id!=''], mp.id[mp.gn!=''&mp.id!='']), unique);
  gn2id <- lapply(split(mp.id[mp.gn!=''&mp.id!=''], mp.gn[mp.gn!=''&mp.id!='']), unique);
  id2gn <- lapply(id2gn, function(x) x[x!='']);
  gn2id <- lapply(gn2id, function(x) x[x!='']);
  id2gn1 <- id2gn[sapply(id2gn, length)>0]; 
  gn2id1 <- gn2id[sapply(gn2id, length)>0]; 
  saveRDS(mp, file=paste(path, 'r/genemap.rds', sep='/'));
  saveRDS(id2gn1, file=paste(path, 'r/genemap_omim2gene.rds', sep='/'));
  saveRDS(gn2id1, file=paste(path, 'r/genemap_gene2omim.rds', sep='/'));

  #######################################################################################
  # genemap2
  mp <- scan(fn.loc[5], sep='\n', flush=TRUE, what='');
  hd <- strsplit(sub('^#', '', mp[4]), '\t')[[1]]; 
  hd <- sub('^ ', '', hd); 
  hd <- sub(' $', '', hd); 
  mp <- mp[substr(mp, 1, 1)!='#'];
  mp <- strsplit(mp, '\t');
  mp <- sapply(mp, function(x) x[1:length(hd)]); 
  mp <- t(mp); 
  colnames(mp) <- hd; 
  mp.gn<-mp[, 'Entrez Gene ID'];
  mp.id<-mp[, 'Mim Number'];
  id2gn <- lapply(split(mp.gn[mp.gn!=''&mp.id!=''], mp.id[mp.gn!=''&mp.id!='']), unique);
  gn2id <- lapply(split(mp.id[mp.gn!=''&mp.id!=''], mp.gn[mp.gn!=''&mp.id!='']), unique);
  id2gn <- lapply(id2gn, function(x) x[x!='']);
  gn2id <- lapply(gn2id, function(x) x[x!='']);
  id2gn2 <- id2gn[sapply(id2gn, length)>0]; 
  gn2id2 <- gn2id[sapply(gn2id, length)>0]; 
  saveRDS(mp, file=paste(path, 'r/genemap2.rds', sep='/'));
  saveRDS(id2gn2, file=paste(path, 'r/genemap2_omim2gene.rds', sep='/'));
  saveRDS(gn2id2, file=paste(path, 'r/genemap2_gene2omim.rds', sep='/'));
  
  #######################################################################################
  # Phenotype series
  tbl <- ImportTable(fn.loc[6]);
  srs <- rownames(tbl);
  names(srs) <- as.vector(tbl[[1]]);
  mps <- lapply(names(srs), function(id) {
    fn <- paste(path, 'src', paste(id, '.html', sep=''), sep='/');
    lk <- paste('https://www.omim.org/phenotypicSeries', id, sep='/');
    download.file(lk, fn);
    tb <- ImportTable(fn, rownames = FALSE, colnames = TRUE);
    colnames(tb) <- tb[1, ]; 
    tb[tb[, 5] %in% rownames(an), , drop=FALSE];
  });
  names(mps) <- names(srs);
  ss2gn <- lapply(mps, function(x) unique(unlist(hu.gn[x[[6]]])));
  ss <- data.frame(Name=srs, Num_Phenotype=sapply(mps, nrow),
                   URL=paste('https://www.omim.org/phenotypicSeries', names(srs), sep='/'));
  rownames(ss) <- names(srs);
  saveRDS(ss, file=paste(path, 'r/phenotype_series.rds', sep='/'));
  saveRDS(mps, file=paste(path, 'r/phenotype_series_full.rds', sep='/'));
  saveRDS(ss2gn, file=paste(path, 'r/phenotype_series2gene.rds', sep='/'));
  
  list(omim=an, series=ss);
}