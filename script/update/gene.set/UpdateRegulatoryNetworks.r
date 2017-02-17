require(rchive);

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/regulatorynetworks', sep='/');
dir.create(path);
dir.create(paste(path, 'r', sep='/'));
dir.create(paste(path, 'src', sep='/'));

library(RCurl);
download.file(
  url = "http://www.regulatorynetworks.org/results/networks/hg19/networks.v09162013.tgz",
  destfile = paste(path, 'src', "networks.v09162013.tgz", sep='/')
);
download.file(
  url = "http://www.regulatorynetworks.org/results/networks/mm9/networks.v12032013.tgz",
  destfile = paste(path, 'src', "networks.v12032013.tgz", sep='/')
);
p1 <- paste(path, 'src', 'human', sep='/');
p2 <- paste(path, 'src', 'mouse', sep='/')
untar(paste(path, 'src', "networks.v09162013.tgz", sep='/'), exdir=p1);
untar(paste(path, 'src', "networks.v12032013.tgz", sep='/'), exdir=p2);

f1 <- dir(p1, rec=TRUE);
f2 <- dir(p2, rec=TRUE);

## Human
library("org.Hs.eg.db", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3");
ks <- mappedkeys(org.Hs.egSYMBOL);
mp <- as.list(org.Hs.egSYMBOL[ks]);
nm <- unlist(mp);
mp <- names(mp); 
names(mp) <- nm; 
f1 <- f1[grep('genes.regulate.genes$', f1)];
c1 <- sapply(strsplit(f1, '/'), function(x) rev(x)[2]);
c1 <- sapply(strsplit(c1, '-'), function(x) x[1]);
m1 <- lapply(paste(p1, f1, sep='/'), function(f) {
  tb <- read.table(f, sep='\t', stringsAsFactors = FALSE);
  id <- mp[tb[, 1]];
  split(as.vector(id[!is.na(id)]), tb[!is.na(id), 2]);
}); 
names(m1) <- c1;
fn <- sapply(c1, function(c) saveRDS(m1[[c]], paste(path, '/r/human_', c, '.rds', sep='')));

df <- data.frame(Cell=rep(names(m1), sapply(m1, length)), TF=unlist(lapply(m1, names)), stringsAsFactors = FALSE);
df$Name <- paste(df$TF, 'targets in human', df$Cell, 'cell');
df$URL <- paste('https://www.ncbi.nlm.nih.gov/gene/?term=', df$TF, sep='');
hs <- do.call('c', m1);
names(hs) <- rownames(df) <- paste(df$Cell, df$TF, sep='-');
saveRDS(df, paste(path, 'r', 'anno_human.rds', sep='/'));
saveRDS(hs, paste(path, 'r', 'mapping_human.rds', sep='/'));

out <- list(human=df); 

## Mouse
library("org.Mm.eg.db", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3");
ks <- mappedkeys(org.Mm.egSYMBOL);
mp <- as.list(org.Mm.egSYMBOL[ks]);
nm <- unlist(mp);
mp <- names(mp); 
names(mp) <- tolower(nm); 
f2 <- f2[grep('genes.regulate.genes$', f2)];
c2 <- sapply(strsplit(f2, '/'), function(x) rev(x)[2]);
c2 <- sapply(strsplit(c2, '-'), function(x) x[1]);
m2 <- lapply(paste(p2, f2, sep='/'), function(f) {
  tb <- read.table(f, sep='\t', stringsAsFactors = FALSE);
  id <- mp[tolower(tb[, 4])];
  split(as.vector(id[!is.na(id)]), tb[!is.na(id), 5]);
}); 
names(m2) <- c2;
m2 <- lapply(unique(c2), function(nm) {
  x <- m2[c2==nm];
  if (length(x)==1) x[[1]] else {
    x <- do.call('c', x); 
    names(x) <- sapply(strsplit(names(x), '\\.'), function(a) a[2]); 
    x <- split(x, names(x));
    lapply(x, function(x) sort(unique(unlist(x, use.names=FALSE))));
  }
});
names(m2) <- unique(c2);
fn <- sapply(c2, function(c) saveRDS(m2[[c]], paste(path, '/r/mouse_', c, '.rds', sep='')));

df <- data.frame(Cell=rep(names(m2), sapply(m2, length)), TF=unlist(lapply(m2, names)), stringsAsFactors = FALSE);
df$Name <- paste("Bound by", df$TF, 'in mouse', df$Cell, 'cell');
df$URL <- paste('https://www.ncbi.nlm.nih.gov/gene/?term=', df$TF, sep='');
mm <- do.call('c', m2);
names(mm) <- rownames(df) <- paste(df$Cell, df$TF, sep='-');
saveRDS(df, paste(path, 'r', 'anno_mouse.rds', sep='/'));
saveRDS(mm, paste(path, 'r', 'mapping_mouse.rds', sep='/'));

out$mouse <- df;

##############################################################################################################
UpdateLog(out, paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/regulatorynetworks', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/UpdateRegulatoryNetworks.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gene.set/log/', tm, '_UpdateRegulatoryNetworks.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)

