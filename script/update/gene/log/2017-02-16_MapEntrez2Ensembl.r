library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

require(biomaRt);

path <- paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r', sep='/');

org <- c(
  'human' = 'org.Hs.eg.db',
  'mouse' = 'org.Mm.eg.db',
  'rat' = 'org.Rn.eg.db',
  'chimp' = 'org.Pt.eg.db',
  'pig' = 'org.Ss.eg.db',
  'chicken' = 'org.Gg.eg.db',
  'dog' = 'org.Cf.eg.db',
  'cow' = 'org.Bt.eg.db',
  'worm' = 'org.Ce.eg.db',
  'fly' = 'org.Dm.eg.db',
  'zebrafish' = 'org.Dr.eg.db',
  'yeast' = 'org.Sc.sgd.db'
); 

mp <- lapply(names(org), function(nm) {
  cat(nm, '\n'); 
  BiocInstaller::biocLite(org[nm], suppressAutoUpdate=TRUE, suppressUpdates=TRUE);
  require(org[nm], character.only=TRUE); 
  en <- sub('.db$', 'ENSEMBL', org[nm]);
  if (exists(en)) {
    x <- get(en); 
    mapped_genes <- mappedkeys(x);
    xx <- as.list(x[mapped_genes]);
    x <- rep(names(xx), sapply(xx, length));
    y <- unlist(xx, use.names = FALSE); 
    cbind(y, x); 
  } else matrix(nr=0, nc=2);
});
names(mp) <- names(org); 

# Mapping info from biomaRt
ds <- c(
  'human' = 'hsapiens_gene_ensembl',
  'mouse' = 'mmusculus_gene_ensembl',
  'rat' = 'rnorvegicus_gene_ensembl',
  'chimp' = 'ptroglodytes_gene_ensembl',
  'pig' = 'sscrofa_gene_ensembl',
  'chicken' = 'ggallus_gene_ensembl',
  'dog' = 'cfamiliaris_gene_ensembl',
  'cow' = 'btaurus_gene_ensembl',
  'worm' = 'celegans_gene_ensembl',
  'fly' = 'dmelanogaster_gene_ensembl',
  'zebrafish' = 'drerio_gene_ensembl',
  'yeast' = 'scerevisiae_gene_ensembl'
);  

mt <- useMart("ensembl");
mp <- lapply(names(ds), function(nm) {
  cat(nm, '\n'); 
  mt <- useDataset(ds[nm], mt); 
  mc <- getBM(attributes = c('ensembl_gene_id', 'entrezgene'), mart=mt);
  mc <- mc[!is.na(mc[, 2]), ]; 
  cbind(c(mc[,1], mp[[nm]][,1]), c(mc[,2], mp[[nm]][, 2]));
});
names(mp) <- names(ds); 

ids <- lapply(names(mp), function(nm) {
  cat(nm, '\n'); 
  mp <- mp[[nm]];
  m1 <- split(mp[, 1], mp[, 2]); 
  m2 <- split(mp[, 2], mp[, 1]); 
  m1 <- lapply(m1, unique);
  m2 <- lapply(m2, unique); 
  f1 <- paste(nm, '_id2ensg.rds', sep='');
  f2 <- paste(nm, '_ensg2id.rds', sep='');
  saveRDS(m1, paste(path, f1, sep='/'));
  saveRDS(m2, paste(path, f2, sep='/'));
  list(Entrez=names(m1), Ensembl=names(m2)); 
});

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/MapEntrez2Ensembl.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/log/', tm, '_MapEntrez2Ensembl.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)

