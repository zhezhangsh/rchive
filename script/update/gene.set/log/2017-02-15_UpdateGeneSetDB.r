require(rchive); 

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesetdb', sep='/');
dir.create(path, showWarnings = FALSE); 
dir.create(paste(path, 'r', sep='/'), showWarnings = FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings = FALSE);

library(RCurl)
ln1 <- readLines('http://www.genesetdb.auckland.ac.nz/download.php?filename=download/gmt_h');
ln2 <- readLines('http://www.genesetdb.auckland.ac.nz/download.php?filename=download/gmt_m');
ln3 <- readLines('http://www.genesetdb.auckland.ac.nz/download.php?filename=download/gmt_r'); 
lns <- list(human=ln1, mouse=ln2, rat=ln3); 

htm <- readLines('http://www.genesetdb.auckland.ac.nz/sourcedb.html');
lnk <- htm[grep('<dt><li>', htm)]; 
lnk <- strsplit(lnk, '"'); 
nms <- CleanHtmlTags(sapply(lnk, function(x) x[5]));
lnk <- sapply(lnk, function(x) x[2]); 
nms <- sapply(strsplit(nms, ' \\('), function(x) x[1]); 
names(nms) <- lnk; 
nms <- nms[nms!="Rel/NF-kappaB target genes"]
nms <- sub(' ', '_', nms); 
nms[nms=='Wikipathways'] <- 'WikiPathways';
nms[nms=='HumnCyc'] <- 'HumanCyc';
lnk <- names(nms); 
names(lnk) <- nms; 

mps <- lapply(lns, function(l) {
  l  <- strsplit(l, '\t');
  nm <- sapply(l, function(l) l[1]); 
  cl <- sapply(l, function(l) l[2]); 
  gn <- sapply(l, function(l) l[3:length(l)]); 
  cl[cl=='KEGG(Disease)'] <- 'KEGG_Disease'
  cl <- gsub(' ', '_', cl); 
  names(gn) <- nm; 
  mp <- split(gn, cl); 
  mp[sapply(mp, length)>1]; 
}); 
url <- lnk[names(mps[[1]])]; 
names(url) <- names(mps[[1]]); 
url[is.na(url)] <- 'http://www.geneontology.org/'; 

saveRDS(url, paste(path, 'r', 'link.rds', sep='/')); 
saveRDS(mps[[1]], paste(path, 'r', 'human', sep='/')); 
saveRDS(mps[[2]], paste(path, 'r', 'mouse', sep='/')); 
saveRDS(mps[[3]], paste(path, 'r', 'rat', sep='/')); 

##############################################################################################################
UpdateLog(out, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/genesetdb', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateGeneSetDB.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateGeneSetDB.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)