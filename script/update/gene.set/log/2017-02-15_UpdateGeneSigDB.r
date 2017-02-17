require(rchive); 

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesigdb', sep='/');
dir.create(path, showWarnings = FALSE); 
dir.create(paste(path, 'r', sep='/'), showWarnings = FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings = FALSE);

sym <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r/human_genes_synonyms2id.rds', sep='/'));

library(RCurl)
lns <- readLines('http://www.genesigdb.org/genesigdb/download/ALL_SIGSv4.gmt');
lns <- strsplit(lns, '\t'); 
lst <- lapply(lns, function(x) x[3:length(x)]); 
lst <- lapply(lst, function(x) unique(as.vector(unlist(sym[x]))));
lst <- lapply(lst, function(x) x[!is.na(x)]); 
pid <- sapply(lns, function(x) x[1]); 
nms <- sapply(lns, function(x) x[2]); 
url <- sapply(strsplit(pid, '-'), function(x) x[1])
url <- paste('https://www.ncbi.nlm.nih.gov/pubmed/?term=', url, sep='');
pid <- paste('PMID', pid, sep='');
ann <- data.frame(Name=nms, URL=url, stringsAsFactors = FALSE); 
rownames(ann) <- names(lst) <- pid; 

saveRDS(lst, paste(path, 'r', 'mapping.rds', sep='/')); 
saveRDS(ann, paste(path, 'r', 'anno.rds', sep='/')); 

##############################################################################################################
UpdateLog(out, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesigdb', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateGeneSigDB.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateGeneSigDB.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE);

