library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

library(yaml);

fn.phred<-yaml.load_file(paste(Sys.getenv('RCHIVE_HOME'), "source/script/update/gwas/phred_tables.yaml", sep='/'));
gwas<-PrepareMetaGwas(unlist(fn.phred), patch.only=TRUE);

# Update log;
UpdateLog(gwas, paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas', sep='/'));

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), 'source/script/update/gwas/UpdateMetaGwas.r', sep='/');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gwas/log/', tm, '_UpdateMetaGwas.r' , sep='');
file.copy(fn0, fn1)