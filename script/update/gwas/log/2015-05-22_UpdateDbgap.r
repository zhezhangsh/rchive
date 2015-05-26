library(devtools);
install_github("zhezhangsh/rchive");
library(rchive);

options(stringsAsFactors=FALSE);

cat('Downloading dbGaP analyses\n');
meta<-DownloadDbGap();
cat('Retrieve dbGaP p values\n');
RetrieveDbGapStat(rownames(meta), meta[,'study'], stat.name='p value');
cat('Summarize dbGaP metadata\n');
SummarizeDbGap(meta);
cat('Add PubMed\n');
AddDbGapPubMed();

cat('Update log\n');
UpdateLog(meta, paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'));

tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), 'source/script/update/gwas/UpdateDbgap.r', sep='/');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/gwas/log/', tm, '_UpdateDbgap.r' , sep='');
file.copy(fn0, fn1)
