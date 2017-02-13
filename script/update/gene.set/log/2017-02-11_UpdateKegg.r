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
pth <- do.call('rbind', pth);
fn0 <- paste(path, '/r/anno_pathway.rds', sep='');
fn1 <- paste(path, '/r/anno_pathway_species.rds', sep='');
saveRDS(pth, );

##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateKegg.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateKegg.r' , sep='');
file.copy(fn0, fn1)
