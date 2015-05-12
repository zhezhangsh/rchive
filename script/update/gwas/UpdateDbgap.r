library(rchive);

meta<-DownloadDbGap();

ana<-RetrieveDbGapStat(rownames(meta), meta[,'study'], stat.name='p value');

UpdateLog(meta, paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'));