library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r', sep='/');

# Selected species
species<-c(
  'human' = '9606',
  'mouse' = '10090',
  'rat' = '10116',
  'chimp' = '9598',
  'pig' = '9823',
  'chicken' = '9031',
  'dog' = '9615',
  'cow' = '9913',
  'worm' = '6239',
  'fly' = '7227',
  'zebrafish' = '7955',
  'ecoli' = '511145',
  'yeast' = '559292'
);
sp <- names(species);
fn <- paste(path, '/', sp, '_genes_full.rds', sep='');
fn <- fn[file.exists(fn)];
id <- lapply(fn, function(f) {
  id <- readRDS(f);
  saveRDS(rownames(id), sub('genes_full', 'just_id_entrez', f));
  saveRDS(as.vector(id$Symbol), sub('genes_full', 'just_id_symbol', f));
});

require("biomaRt");
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
ids <- lapply(names(ds), function(nm) {
  cat(nm, '\n'); 
  mt <- useDataset(ds[nm], mt); 
  mc <- getBM(attributes = c('ensembl_gene_id'), mart=mt);
  id <- unique(mc[, 1]); 
  id <- id[!is.na(id)];
  fn <- paste(path, '/', nm, '_just_id_ensembl.rds', sep='');
  saveRDS(id, fn); 
  id;
});

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/UpdateGeneID.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/log/', tm, '_UpdateGeneID.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)
