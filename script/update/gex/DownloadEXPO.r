devtools::install_github("zhezhangsh/rchive");

library(rchive);
library(affy);
library(devtools);

download.all<-FALSE;

options(stringsAsFactors=FALSE);

path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/expo', sep='/');

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

# path.geodb<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/db', sep='/');
# if (!file.exists(path.geodb)) dir.create(path.geodb);
# fn.geodb<-paste(path.geodb, 'GEOmetadb.sqlite', sep='/');
# if (!file.exists(fn.geodb) | download.all) GEOmetadb::getSQLiteFile(path.geodb);
# 
# # Get sample list
# gse.id<-'GSE2109';
# db<-dplyr::src_sqlite(fn.geodb);
# mp<-dplyr::tbl(db, 'gse_gsm');
# gsm.id<-as.data.frame(dplyr::collect(dplyr::filter(mp, gse==gse.id)))[, 1];

# download series metadata
fn.meta<-paste(path, 'r', 'original_metadata.rds', sep='/');
if (!file.exists(fn.meta) | download.all) {
  gse<-getGEO(gse.id, GSEMatrix=FALSE, AnnotGPL=FALSE, getGPL=FALSE);
  meta<-list(GSE=gse@header, GSMs=lapply(gse@gsms, function(g) g@header)); 
  saveRDS(meta, file=fn.meta);
} else meta<-readRDS(fn.meta);

# sample metadata
gsms<-meta[[2]];
desc<-lapply(gsms, function(g) {
  d<-g$description;
  d<-strsplit(d, ':');
  v<-sapply(d, function(d) d[2]);
  v<-sub('^ ', '', v);
  v<-sub('\\s+', '', v);
  names(v)<-sapply(d, function(d) d[1]);
  v;
});
cnms<-sort(unique(unlist(lapply(desc, names), use.names=FALSE)));


# download CEL files
gsm.id<-names(meta$GSMs);
fn.cel<-paste(path, '/src/', gsm.id, '.CEL.gz', sep='');
names(fn.cel)<-gsm.id;
sapply(gsm.id, function(id) if (download.all | !file.exists(fn.cel[id])) 
  GEOquery::getGEOSuppFiles(id, FALSE, paste(path, 'src', sep='/')) else NA)->na;



