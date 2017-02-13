library(devtools);
# source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");
install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

# download source files from PubTator FTP server
ftp.files<-c(
  chemical="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//chemical2pubtator.gz",
  species="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//species2pubtator.gz",
  disease="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//disease2pubtator.gz",
  gene="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//gene2pubtator.gz",
  mutation="ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator//mutation2pubtator.gz"
)
for (i in 1:length(ftp.files)) {
  cat('Loading', ftp.files[i]);
  ParsePubtator(ftp.files[i]);
}

fn  <- paste(path, '/r/original_pubmed2gene.rds', sep='')
p2g <- readRDS(fn);
pid <- p2g[,1];
gid <- p2g[,2];
gid0<- strsplit(gid, ',');
n   <-sapply(gid0, length);
pid <- rep(pid, n);
gid <- unlist(gid0, use.names=FALSE);
mp1 <- split(gid, pid);
mp2 <- split(pid, gid); 
mp1 <- lapply(mp1, unique);
mp2 <- lapply(mp2, unique);

path = paste(Sys.getenv("RCHIVE_HOME"), "data/literature/public/pubtator/r", sep = "/")
saveRDS(mp1, paste(path, 'map_pubmed2gene.rds', sep='/'));
saveRDS(mp2, paste(path, 'map_gene2pubmed.rds', sep='/'));

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), 'source/script/update/literature/UpdatePubtator.r', sep='/');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/literature/log/', tm, '_UpdatePubtator.r' , sep='');
file.copy(fn0, fn1)

