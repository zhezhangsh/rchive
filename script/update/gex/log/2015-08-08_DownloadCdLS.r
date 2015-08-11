devtools::install_github("zhezhangsh/rchive");

library(rchive);

options(stringsAsFactors=FALSE);

path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/cdls', sep='/');

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

geo<-c("GSE55316", "GSE64034", "GSE12408", "GSE21844");
paths<-sapply(geo, function(id) ProcessGEO(id, paste(path, 'src', sep='/'))); 
names(paths)<-geo;

##############################################################################################################
UpdateLog(paths, paste(Sys.getenv("RCHIVE_HOME"), 'data/gex/public/', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/DownloadCdLS.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gex/log/', tm, '_DownloadCdLS.r' , sep='');
file.copy(fn0, fn1)
