library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

# Types of KEGG databases
types <- c('pathway', 'module', 'compound', 'reaction', 'enzyme', 'disease', 'drug', 'dgroup', 'omim', 'pfam');

# Species to map
species <- c('human'='hsa', 'mouse'='mmu', 'rat'='rno', 'fly'='dme', 'worm'='cel', 'zebrafish'='dre', 'yeast'='sce', 'ecoli'='eco');

# Path to output files
out <- paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg', sep='/'); 

# Download data and save results
ids <- MapKegg2Gene.2(types, species, out); 

# Pathway by species;
url <- paste('http://rest.kegg.jp/list/pathway', species, sep='/');
names(url) <- names(species);
url <- url[sapply(url, url.exists)];
pth <- lapply(names(url), function(nm) { 
  cat(nm, '\n'); 
  l <- readLines(url[nm]); 
  l <- strsplit(l, '\t');
  x <- sapply(l, function(x) x[1]);
  y <- sapply(l, function(x) x[2]);
  x <- sapply(strsplit(x, ':'), function(x) x[2]); 
  names(y) <- x;
  df <- data.frame(Species=rep(nm, length(x)), Name=y, 
                   URL=paste('http://www.kegg.jp/dbget-bin/www_bget?', x, sep=''), stringsAsFactors = FALSE);
  rownames(df) <- x;
  df; 
});
fn0 <- paste(out, '/r/anno_pathway.rds', sep='');
fn1 <- paste(out, '/r/anno_pathway_common.rds', sep='');
file.rename(fn0, fn1);
pth <- do.call('rbind', pth);
pth <- pth[, -1];
pth[, 1] <- sapply(strsplit(pth[, 1], ' - '), function(x) x[1]);
saveRDS(pth, fn0);

# Module by species;
url <- paste('http://rest.kegg.jp/list/module', species, sep='/');
names(url) <- names(species);
url <- url[sapply(url, url.exists)];
pth <- lapply(names(url), function(nm) { 
  cat(nm, '\n'); 
  l <- readLines(url[nm]); 
  l <- strsplit(l, '\t');
  x <- sapply(l, function(x) x[1]);
  y <- sapply(l, function(x) x[2]);
  x <- sapply(strsplit(x, ':'), function(x) x[2]); 
  names(y) <- x;
  df <- data.frame(Species=rep(nm, length(x)), Name=y, 
                   URL=paste('http://www.kegg.jp/dbget-bin/www_bget?', x, sep=''), stringsAsFactors = FALSE);
  rownames(df) <- x;
  df; 
});
fn0 <- paste(out, '/r/anno_module.rds', sep='');
fn1 <- paste(out, '/r/anno_module_common.rds', sep='');
file.rename(fn0, fn1);
pth <- do.call('rbind', pth);
pth <- pth[, -1];
saveRDS(pth, fn0);

# PFAM
require(PFAM.db);
x <- PFAMDE;
mapped_keys <- mappedkeys(x); 
xx <- as.list(x[mapped_keys]); 
id <- names(xx); 
nm <- unlist(xx, use.names=FALSE); 
df <- data.frame(Name=nm, URL=paste("http://pfam.xfam.org/family", id, sep='/')); 
rownames(df) <- id;
saveRDS(df, paste(out, '/r/anno_pfam.rds', sep=''));

# OMIM
mm <- readRDS('/zhangz/rchive/data/disease/public/omim/r/omim_anno.rds')
mm <- mm[, c('Preferred_Title', 'URL')];
names(mm)[1] <- 'Name';
saveRDS(mm, paste(out, '/r/anno_omim.rds', sep=''));

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateKegg.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateKegg.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)
