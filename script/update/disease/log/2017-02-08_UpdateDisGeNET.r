library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

lnk <- c(
  curated = 'http://www.disgenet.org/ds/DisGeNET/results/curated_gene_disease_associations.tsv.gz',
  literature = "http://www.disgenet.org/ds/DisGeNET/results/literature_gene_disease_associations.tsv.gz", 
  befree = "http://www.disgenet.org/ds/DisGeNET/results/befree_gene_disease_associations.tsv.gz",
  all = "http://www.disgenet.org/ds/DisGeNET/results/all_gene_disease_associations.tsv.gz"
); 

dis<-ParseDisGeNET(url = lnk);

##############################################################################################################
UpdateLog(dis, paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/disgenet', sep='/'), just.new=FALSE);

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/disease/UpdateDisGeNET.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/disease/log/', tm, '_UpdateDisGeNET.r' , sep='');
file.copy(fn0, fn1)
