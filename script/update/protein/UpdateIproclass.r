library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

########################################################################################
download.again <- FALSE;
path <- paste(Sys.getenv("RCHIVE_HOME"), "data/protein/public/iproclass", sep = "/");
########################################################################################

tax <- c(
  'human' = '9606',
  'mouse' = '10090',
  'rat' = '10116',
  'chimp' = '9598',
  'pig' = '9823',
  'chicken' = '9031',
  'dog' = '9615',
  'cow' = '9913',
  'worm' = '6239',
  'fly' = '7227',
  'zebrafish' = '7955',
  'yeast' = '559292'
);

if (download.again) {
  setwd(paste(path, 'src', sep='/'));
  system(paste('sh', paste(path, 'src', 'shell.sh', sep='/')));
}; 

tx <- readLines(paste(path, 'src', 'column_tax.tsv', sep='/')); 
ind <- which(tx %in% tax);
gn <- readLines(paste(path, 'src', 'column_gene.tsv', sep='/'))[ind]; 
pf <- readLines(paste(path, 'src', 'column_pfam.tsv', sep='/'))[ind]; 
go <- readLines(paste(path, 'src', 'column_go.tsv', sep='/'))[ind]; 
sf <- readLines(paste(path, 'src', 'column_pirsf.tsv', sep='/'))[ind]; 
om <- readLines(paste(path, 'src', 'column_omim.tsv', sep='/'))[ind]; 
pm <- readLines(paste(path, 'src', 'column_pubmed.tsv', sep='/'))[ind]; 
tx <- tx[ind];

all <- list(pfam=pf, go=go, pirsf=sf, omim=om, pubmed=pm); 
ids <- lapply(names(all), function(nm1) {
  cat(nm1, '\n'); 
  x <- all[[nm1]]; 
  lapply(names(tax), function(nm2) {
    cat(nm2, '\n');
    i <- which(tx==tax[nm2] & x!='' & gn!='');
    y <- x[i];
    g <- gn[i];
    if (length(g) > 0) {
      z <- strsplit(y, '; '); 
      y <- rep(g, sapply(z, length)); 
      z <- unlist(z, use.names=FALSE); 
      y <- strsplit(y, '; ');
      z <- rep(z, sapply(y, length));
      y <- unlist(y, use.names=FALSE); 
      z <- sub('/$', '', z);
      m <- split(y, z); 
      m <- lapply(m, unique);
      m <- m[sapply(m, length)>0];
      f <- paste(path, 'r', paste(nm2, '_', nm1, '2gene.rds', sep=''), sep='/');
      saveRDS(m, file=f); 
      names(m); 
    } else c();
  });
});

ids <- lapply(ids, unlist); 
ids <- lapply(ids, unique);
ids <- lapply(ids, sort); 
names(ids) <- names(all);

##############################################################################################################
# Annotation
## PFAM
fn1 <- paste(Sys.getenv("RCHIVE_HOME"), "data/protein/public/iproclass/src/Pfam-A.clans.tsv.gz", sep = "/");
pfm <- read.csv(fn1, sep='\t', header = FALSE, row.names = 1, stringsAsFactors = FALSE);
colnames(pfm) <- c('Clan', 'Clan_Symbol', 'Symbol', 'Name');
pfm <- pfm[, c(3, 4, 1, 2)];
pfm$URL <- paste('http://pfam.xfam.org/family', rownames(pfm), sep='/');
saveRDS(pfm, file=paste(path, 'r', 'anno_pfam.rds', sep='/')); 
# pfam clan
fn2 <- paste(Sys.getenv("RCHIVE_HOME"), "data/protein/public/iproclass/src/Pfam-C.gz", sep = "/");
lns <- readLines(fn2);
acc <- lns[grep('^#=GF AC', lns)];
acc <- strsplit(acc, '[ \\.]');
acc <- sapply(acc, function(x) rev(x)[2]);
nam <- lns[grep('^#=GF DE', lns)];
nam <- sub('^#=GF DE   ', '', nam);
cln <- data.frame(Name=nam, URL=paste('http://pfam.xfam.org/clan', acc, sep='/'), stringsAsFactors = FALSE);
rownames(cln) <- acc;
saveRDS(cln, file=paste(path, 'r', 'anno_pfamclan.rds', sep='/')); 
# clan2gene
fns <- dir(paste(path, 'r', sep='/'));
fns <- fns[grep('pfam2gene', fns)];
fns <- paste(path, 'r', fns, sep='/');
id0 <- lapply(fns, function(f) {
  m <- readRDS(f); 
  p <- pfm[pfm$Clan %in% rownames(cln), , drop=FALSE];
  m <- m[names(m) %in% rownames(p)];
  m <- split(m, p[names(m), 'Clan']); 
  m <- lapply(m, function(m) unique(unlist(m, use.names=FALSE))); 
  saveRDS(m, file=sub('pfam2gene', 'pfamclan2gene', f)); 
  names(m); 
});
ids$pfamclan <- unique(unlist(id0));

## GO
library(GO.db);
fns <- dir(paste(path, 'r', sep='/'));
fns <- fns[grep('go2gene', fns)];
fns <- paste(path, 'r', fns, sep='/');
id0 <- lapply(fns, function(f) names(readRDS(f))); 
id0 <- sort(unique(unlist(id0, use.names=FALSE)));
got <- as.list(GOTERM);
got <- got[names(got) %in% id0];
nam <- sapply(got, function(x) x@Term);
ann <- data.frame(Name=nam, URL=paste('http://amigo.geneontology.org/amigo/term', names(got), sep='/'), 
                  stringsAsFactors = FALSE);
rownames(ann) <- names(got);
saveRDS(ann, file=paste(path, 'r', 'anno_go.rds', sep='/')); 

## PIRSF
fns <- dir(paste(path, 'r', sep='/'));
fns <- fns[grep('pirsf2gene', fns)];
fns <- paste(path, 'r', fns, sep='/');
id0 <- lapply(fns, function(f) names(readRDS(f))); 
id0 <- sort(unique(unlist(id0, use.names=FALSE)));
fn1 <- paste(Sys.getenv("RCHIVE_HOME"), "data/protein/public/iproclass/src/pirsfinfo.dat", sep = "/");
lns <- readLines(fn1);
lns <- lns[grep('^>', lns)];
id0 <- sub('^>', '', sapply(strsplit(lns, ' '), function(x) x[1]));
nam <- sapply(strsplit(lns, '\\) '), function(x) x[2]);
nam[is.na(nam)] <- '';
ann <- data.frame(Name=nam, URL=paste('http://pir.georgetown.edu/cgi-bin/ipcSF?id=', id0, sep=''), 
                  stringsAsFactors = FALSE);
rownames(ann) <- id0; 
saveRDS(ann, file=paste(path, 'r', 'anno_pirsf.rds', sep='/')); 

## OMIM
omm <- readRDS(paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/omim/r/omim_anno.rds', sep='/'));
id0 <- names(readRDS(paste(path, 'r', 'human_omim2gene.rds', sep='/')));
ann <- omm[id0, c('Preferred_Title', 'URL')];
names(ann)[1] <- 'Name';
saveRDS(ann, file=paste(path, 'r', 'anno_omim.rds', sep='/')); 

## PubMed
fns <- dir(paste(path, 'r', sep='/'));
fns <- fns[grep('pubmed2gene', fns)];
fns <- paste(path, 'r', fns, sep='/');


##############################################################################################################
UpdateLog(ids, paste(Sys.getenv("RCHIVE_HOME"), 'data/protein/public/iproclass', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/disease/UpdateIproclass.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/disease/log/', tm, '_UpdateIproclass.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE)
