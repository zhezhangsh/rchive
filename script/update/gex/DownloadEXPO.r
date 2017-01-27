devtools::install_github("zhezhangsh/rchive");

library(rchive);
library(affy);
library(devtools);

download.all<-FALSE;

options(stringsAsFactors=FALSE);

path   <- paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/expo', sep='/');
gse.id <- 'GSE2109';
cdf.nm <- 'hgu133plus2hsentrezg'; 

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

# download series metadata
fn.meta<-paste(path, 'r', 'original_metadata.rds', sep='/');
if (!file.exists(fn.meta) | download.all) {
  gse <- GEOquery::getGEO(gse.id, GSEMatrix=FALSE, AnnotGPL=FALSE, getGPL=FALSE);
  meta<-list(GSE=gse@header, GSMs=lapply(gse@gsms, function(g) g@header)); 
  saveRDS(meta, file=fn.meta);
} else meta<-readRDS(fn.meta);
saveRDS(meta$GSMs, file=paste(path, 'r', 'original_metadata.rds', sep='/'));

# sample metadata
gsms <- meta[[2]];
desc <- lapply(gsms, function(g) {
  d<-g$description;
  d<-strsplit(d, ':');
  v<-sapply(d, function(d) d[2]);
  v<-sub('^ ', '', v);
  v<-sub('\\s+', '', v);
  names(v)<-sapply(d, function(d) d[1]);
  v;
});
cnms <- sort(unique(unlist(lapply(desc, names), use.names=FALSE)));
cnms <- sort(Reduce('union', cnms));
desc <- t(sapply(desc, function(x) x[cnms])); 
desc <- data.frame(desc, stringsAsFactors = FALSE);
dimnames(desc) <- list(names(gsms), cnms); 
ttls <- sapply(gsms, function(x) x$title); 
ttls <- gsub(' - ', '-', ttls);
ttls <- as.vector(gsub(' ', '_', ttls)); 
prim <- as.vector(desc$`Primary Site`); 
prim <- gsub(' ', '_', prim); 
qual <- as.numeric(desc[, 'Quality metric = 28S to 18S']);
gnd  <- desc[, 'Gender'];
age  <- desc[, 'Patient Age'];
toba <- desc[, 'Tobacco Use '];
ethc <- desc[, 'Ethnic Background'];
fami <- desc[, 'Family History of Cancer?'];
alco <- desc[, 'Alcohol Consumption?'];

# download CEL files
gsm.id<-names(meta$GSMs);
fn.cel<-paste(path, '/src/', gsm.id, '.CEL.gz', sep='');
names(fn.cel)<-gsm.id;
for (i in 1:length(fn.cel)) {
  print(i); 
  if (download.all | !file.exists(fn.cel[i]))
    try(GEOquery::getGEOSuppFiles(names(fn.cel)[i], FALSE, paste(path, 'src', sep='/')));
}

# Re-try if failed earlier
fn0 <- fn.cel[!file.exists(fn.cel)]; 
while(length(fn0) > 0) {
  for (i in 1:length(fn0)) {
    print(i); 
    if (download.all | !file.exists(fn0[i]))
      try(GEOquery::getGEOSuppFiles(names(fn0)[i], FALSE, paste(path, 'src', sep='/')));
  }
}

# Load cel file and normalize data
raw <- ReadAffy(filenames = fn.cel);
raw@cdfName <- cdf.nm; 
expr <- exprs(rma(raw)); 
colnames(expr) <- sub('.CEL.gz', '', colnames(expr)); 
rownames(expr) <- sub('_at', '', rownames(expr)); 
saveRDS(expr, paste(path, 'src', 'expr_all.rds', sep='/')); 



