library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

path <- paste(Sys.getenv("RCHIVE_HOME"), 'data/noncoding/public/rnacentral', sep='/');
path0 <- paste(path, 'r', sep='/'); 
path1 <- paste(path, 'src/current_release/genome_coordinates', sep='/');
path2 <- paste(path, 'src/current_release/id_mapping', sep='/');

gnm <- readRDS(paste(path, 'src/genome.rds', sep='/'));
for (i in 1:nrow(gnm)) {
  print(i); 
  fin <- paste(path1, paste(gnm[i, 1], gnm[i, 2], 'gff.gz', sep='.'), sep='/');
  fout <- paste(rownames(gnm)[i], '_', gnm[i, 2], '_gff.rds', sep='');
  rchive::ParseGtf(fin, fout, path=path0, Dbxref = FALSE)
};

mpp <- read.table(paste(path2, 'id_mapping.tsv.gz', sep='/'), sep='\t', header = FALSE, stringsAsFactors = FALSE);
mpp[, 4] <- as.character(mpp[, 4]); 
names(mpp) <- c('ID', 'Source', 'Sequence', 'Species', 'Type', 'Name');
saveRDS(mpp, paste(path0, 'id_mapping.rds', sep='/'));

# hsa <- mpp[mpp$Species=='9606', ];
# saveRDS(hsa, paste(path0, 'id_mapping_human.rds', sep='/'));

tax <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), "data/taxonomy/public/ncbi/r/name_representative.rds", sep='/'));
tax <- tax[names(tax) %in% mpp$Species]; 

tax <- c('9606'='human', '10116'='rat', '10090'='mouse', '9598'='chimp');
spe <- lapply(names(tax), function(tax) mpp[mpp$Species==tax, -4]); 
names(spe) <- tax;
fns <- sapply(names(spe), function(nm) {
  f <- paste(path0, paste('id_mapping_', nm, '.rds', sep=''), sep='/');
  saveRDS(spe[[nm]], f);
  f;
});

##############################################################################################################
# UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/unigene', sep='/'), just.new=FALSE);
# 
# tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
# fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/UpdateUnigene.r', sep='');
# fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene/log/', tm, '_UpdateUnigene.r' , sep='');
# file.copy(fn0, fn1); 
