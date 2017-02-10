library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);
library(NCBI2R);

options(stringsAsFactors=FALSE);

#url<-"https://www.ebi.ac.uk/gwas/api/search/downloads/full";

updates<-ParseGwasCatalog();

##############################################################################################################
UpdateLog(updates, paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/gwascatalog', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), 'source/update/gwas/UpdateGwasCatalog.r', sep='/');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/update/gwas/log/', tm, '_UpdateGwasCatalog.r' , sep='');
file.copy(fn0, fn1)
