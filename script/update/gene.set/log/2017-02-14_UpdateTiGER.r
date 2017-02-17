# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/tiger', sep='/');
dir.create(path, showWarnings = FALSE); 
dir.create(paste(path, 'r', sep='/'), showWarnings = FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings = FALSE);

library(RCurl)
library('org.Hs.eg.db');
x <- org.Hs.egREFSEQ
# Get the entrez gene identifiers that are mapped to any RefSeq ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
x <- rep(names(xx), sapply(xx, length));
names(x) <- unlist(xx, use.names=FALSE);

lns <- readLines('http://bioinfo.wilmer.jhu.edu/tiger/download/ref2tissue-Table.txt')[-1];
lns <- strsplit(lns, '\t'); 
ref <- sapply(lns, function(x) x[1]); 
tis <- lapply(lns, function(x) x[-1]); 

ref <- x[ref]; 
tis <- tis[!is.na(ref)]; 
ref <- as.vector(ref[!is.na(ref)]); 
ref <- rep(ref, sapply(tis, length)); 
tis <- unlist(tis, use.names=FALSE); 
map <- split(ref, tis); 
saveRDS(map, paste(path, 'r', 'tissue2gene.rds', sep='/'));

##############################################################################################################
UpdateLog(out, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/tiger', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateTiGER.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateTiGER.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)