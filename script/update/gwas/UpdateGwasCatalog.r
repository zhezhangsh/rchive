library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

#url<-"https://www.ebi.ac.uk/gwas/api/search/downloads/full";

updates<-ParseGwasCatalog(update.all.pubmed=TRUE);

##############################################################################################################
UpdateLog(updates, paste(RCHIVE_HOME, 'data/gwas/public/gwascatalog', sep='/'), just.new=TRUE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(RCHIVE_HOME, 'source/update/gwas/UpdateGwasCatalog.r', sep='/');
fn1<-paste(RCHIVE_HOME, '/source/update/gwas/log/', tm, '_UpdateGwasCatalog.r' , sep='');
file.copy(fn0, fn1)
