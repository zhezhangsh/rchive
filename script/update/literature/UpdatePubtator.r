library(devtools);
# source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");
install_github("zhezhangsh/rchive");
library(rchive);

path <- paste(Sys.getenv('RCHIVE_HOME'), 'data/literature/public/pubtator', sep='/');

options(stringsAsFactors=FALSE);

# download source files from PubTator FTP server
ftp.files<-c(
  chemical="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//chemical2pubtator.gz",
  species="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//species2pubtator.gz",
  disease="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//disease2pubtator.gz",
  gene="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//gene2pubtator.gz",
  mutation="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//mutation2pubtator.gz"
)
if (FALSE) {
  for (i in 1:length(ftp.files)) {
    cat('Loading', ftp.files[i]);
    ParsePubtator(ftp.files[i]);
  }
}

########################################################################################################
fn  <- paste(path, '/r/gene2pubmed_original.rds', sep='')
p2g <- readRDS(fn);
pid <- p2g[,1];
gid <- p2g[,2];
src <- p2g[,4];
gid0<- strsplit(gid, '[,;]');
n   <-sapply(gid0, length);
pid <- rep(pid, n);
src <- rep(src, n); 
gid <- unlist(gid0, use.names=FALSE);
mp1 <- split(gid, pid);
mp2 <- split(pid, gid); 
mp1 <- lapply(mp1, unique);
mp2 <- lapply(mp2, unique);

saveRDS(mp1, paste(path, 'r', 'pubmed2gene_all.rds', sep='/'));
saveRDS(mp2, paste(path, 'r', 'gene2pubmed_all.rds', sep='/'));

# Split by species and source
path.gene <- paste(Sys.getenv('RCHIVE_HOME'), 'data/gene/public/entrez/r', sep='/');
species <- c("human", "mouse", "rat", "chimp", "pig", "chicken", "dog", "cow", "worm", "fly", "zebrafish", "yeast", "ecoli"); 
fn.gene <- paste(path.gene, '/', species, '_genes_full.rds', sep='');
names(fn.gene) <- species;
fn.gene <- fn.gene[file.exists(fn.gene)];
gns <- lapply(fn.gene, function(f) rownames(readRDS(f)));
names(gns) <- names(fn.gene);
fns <- lapply(names(gns), function(sp) {
  ind <- which(gid %in% gns[[sp]]); 
  
  pd <- pid[ind];
  gd <- gid[ind];
  sr <- src[ind]; 
  m1 <- split(gd, pd);
  m2 <- split(pd, gd); 
  m1 <- lapply(m1, unique);
  m2 <- lapply(m2, unique);
  saveRDS(m1, paste(path, 'r', paste('pubmed2gene_all_', sp, '.rds', sep=''), sep='/'));
  saveRDS(m2, paste(path, 'r', paste('gene2pubmed_all_', sp, '.rds', sep=''), sep='/'));

  fn <- paste(path, 'r', paste('gene'))
  sc <- strsplit(sr, '\\|'); 
  sc <- sort(unique(unlist(sc, use.names=FALSE)));
  sc <- unique(gsub(' ', '', sc)); 
  
  fn <- sapply(sc, function(nm) {
    cat(nm, '\t'); 
    i <- grep(nm, sr);
    p <- pd[i];
    g <- gd[i]; 
    m <- lapply(split(g, p), unique); 
    f <- paste(path, 'r', paste('pubmed2gene_', sp, '_', nm, '.rds', sep=''), sep='/');
    saveRDS(m, f); 
    f;
  });
  
  cat('\n', sp, length(ind), length(sc), '\n'); 
  fn;
});

########################################################################################################################
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), 'source/script/update/literature/UpdatePubtator.r', sep='/');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/literature/log/', tm, '_UpdatePubtator.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)

