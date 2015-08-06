#devtools::install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

path.coll<-c(
  'ADHB'=paste(Sys.getenv('RCHIVE_HOME'), 'data/gex/public/adhb/r', sep='/'),
  'MAGE'=paste(Sys.getenv('RCHIVE_HOME'), 'data/gex/public/mage/r', sep='/'),
  'Demo'=paste(Sys.getenv('RCHIVE_HOME'), 'data/gex/public/demo/r', sep='/'),
  'ToMD'=paste(Sys.getenv('RCHIVE_HOME'), 'data/gex/public/tomd/r', sep='/')
);

mapping<-lapply(path.coll, PrepareGeexCollection);

##############################################################################################################
UpdateLog(mapping, paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/UpdateGeexCollections.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/log/', tm, '_UpdateGeexCollections.r' , sep='');
file.copy(fn0, fn1)
