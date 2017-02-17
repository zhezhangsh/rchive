library(rchive)

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/encode/public/tfbs', sep='/');
dir.create(path);
dir.create(paste(path, 'r', sep='/'));
dir.create(paste(path, 'src', sep='/'));

library(RCurl);
library(GenomicsRanges);
library(TxDb.Hsapiens.UCSC.hg19.knownGene);

download.file(
  url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz",
  destfile = paste(path, 'src', "wgEncodeRegTfbsClusteredWithCellsV3.bed.gz", sep='/')
);
tbl <- read.table(paste(path, 'src', 'wgEncodeRegTfbsClusteredWithCellsV3.bed.gz', sep='/'), stringsAsFactors = FALSE)
rownames(tbl) <- 1:nrow(tbl); 
loc <- GRanges(tbl[, 1], IRanges(tbl[, 2], tbl[, 3]));

txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, columns=c('gene_id'));
tss <- resize(txs, 1);
tss <- resize(tss, 2000, fix='end');
tss <- resize(tss, 2500, fix='start');
tss <- tss[sapply(tss$gene_id, length)>0];

olp <- as.matrix(findOverlaps(loc, tss));
tbl <- tbl[unique(olp[, 1]), ];
tbl <- tbl[tbl[, 5]>300, ];
gns <- as.list(tss$gene_id);
gns <- gns[olp[, 2]];
gns <- split(gns, olp[, 1]);
gns <- lapply(gns, unlist);
gns <- lapply(gns, unique);
gns <- gns[as.character(unique(olp[, 1]))]

cll <- tbl[, 6]; 
cll <- strsplit(cll, ','); 
ids <- rep(rownames(tbl), sapply(cll, length)); 
map <- split(ids, unlist(cll, use.names=FALSE)); 

all <- lapply(names(map), function(nm1) {
  cat(nm1, '\n'); 
  i <- map[[nm1]]; 
  t <- tbl[i, 4];
  g <- gns[i];
  m <- split(g, t);
  m <- lapply(m, function(m) unique(as.vector(unlist(m, use.names=FALSE))));
  f <- paste(path, 'r', paste('tf_targets_', nm1, '.rds', sep=''), sep='/');
  saveRDS(m, f);
  m; 
}); 
names(all) <- names(map); 
ids <- paste(unlist(lapply(all, names)), rep(names(all), sapply(all, length)), sep='_');
nms <- paste(unlist(lapply(all, names)), 'targets in', rep(names(all), sapply(all, length)));
all <- do.call('c', all); 
names(all) <- ids;

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/encode/public/tfbs', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/encode/MapEncodeTFBS.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/encode/log/', tm, '_MapEncodeTFBS.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE);

