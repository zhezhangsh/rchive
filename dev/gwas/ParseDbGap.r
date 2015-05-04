# Download and process GWAS results from dbGap
DownloadDbGap<-function(url="ftp://ftp.ncbi.nlm.nih.gov/dbgap//Analysis_Table_of_Contents.txt", 
                     path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'), 
                     update.all.analysis=FALSE,
                     redundant.studies=c('216', '88', '342')) {
  # url                   File with analysis metadata and download location
  # path                  Path to output files
  # update.all.analysis   Re-download all PubMed entries if TRUE, takes very long time to do so
  # redundant.studies     Exclude redundant studies of the same analyses belong to other studys
  
  library(RCurl);
  library(NCBI2R);
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if (!file.exists(paste(path, 'r/full_table', sep='/'))) dir.create(paste(path, 'r/full_table', sep='/'), recursive=TRUE);
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
    
  # Original metadata with each ftp file of analysis; parsed into a table
  fn.meta<-paste(path, 'r', 'original_metadata.rds', sep='/');
  if (file.exists(fn.meta)) meta<-readRDS(fn.meta) else {
    emp<-character(0);
    meta<-data.frame(name=emp, study=emp, genome=emp, dbsnp=emp, description=emp, method=emp, file=emp, stringsAsFactors=FALSE);
    saveRDS(meta, file=fn.meta);
  }
  # Original metadata with each ftp file of analysis; as a list
  fn.meta.lst<-paste(path, 'r', 'original_metadata_list.rds', sep='/');
  if (file.exists(fn.meta.lst)) meta.lst<-readRDS(fn.meta.lst) else {
    meta.lst<-lst();
     saveRDS(meta.lst, file=fn.meta.lst);
  }
  
  # Download and parse metadata table
  ln<-strsplit(getURL(url), '\n')[[1]];
  gw<-do.call('rbind', strsplit(ln, '\t'));
  colnames(gw)<-gw[1,];
  gw<-gw[-1, , drop=FALSE];
  gw<-data.frame(gw, stringsAsFactors=FALSE);
  gw<-gw[!(gw[,1] %in% redundant.studies), ];
  gw<-gw;
  
  ######################################################################################
  # Download results of individual analyses
  # file information
  fn.ftp<-gw$analysis.result;
  fn.gz<-sapply(strsplit(fn.ftp, '/'), function(x) x[length(x)]);
  id.study<-sapply(strsplit(fn.gz, '\\.'), function(x) x[1]);
  id.analysis<-sapply(strsplit(fn.gz, '\\.'), function(x) x[2]);
  fn.gz<-paste(path, 'src', fn.gz, sep='/');
  fn.r<-paste(path, 'r/full_table', paste(id.analysis, '.rds', sep=''), sep='/');
  file.info<-cbind(analysis=id.analysis, study=id.study, ftp=fn.ftp, gz=fn.gz, r=fn.r);
  
  ####################################################################################
  # load in table and save to local
  if (update.all.analysis) no.tbl<-file.info else no.tbl<-file.info[!file.exists(file.info[, 'r']), , drop=FALSE];
  if (nrow(no.tbl) > 0) {
    hd<-lapply(1:nrow(no.tbl), function(i) {
      x<-no.tbl[i,];
      cat("loading analysis", x[1], '\n');
      if (!file.exists(x[4])) download.file(x[3], x[4]); 
      ln<-scan(x[4], what='', flush=TRUE, sep='\n');
      
      # very tediously parse file header and add header info to metadata table and list
      hd<-ln[grep('^# ', ln)];
      hd<-sub('^# ', '', hd);
      hd<-strsplit(hd, ':\t');
      hd<-lapply(hd, function(hd) sub('^ ', '', hd));
      hd<-lapply(hd, function(hd) sub(' $', '', hd));      
      hd<-do.call('rbind', hd);
      hd.nm<-tolower(hd[,1]);
      hd<-hd[,2];
      names(hd)<-hd.nm;
      
      # save full result table
      cnm<-strsplit(ln[length(hd)+1], '\t')[[1]];
      tbl<-strsplit(ln[(length(hd)+2):length(ln)], '\t');
      tbl<-t(sapply(tbl, function(x) x[1:length(cnm)]));
      tbl<-data.frame(tbl, stringsAsFactors=FALSE); # all columns are character
      colnames(tbl)<-cnm;
      saveRDS(tbl, x[5]);
      
      hd;
    })
    
    meta.lst[no.tbl[,1]]<-hd;
    new.rows<-sapply(1:length(hd), function(i) {
      hd<-hd[[i]];
      as.vector(c(hd['name'], no.tbl[i, 2], hd[c('human genome build', 'dbsnp build', 'description', 'method')], no.tbl[i, 3]));
    });
    new.rows<-t(new.rows);
    rownames(new.rows)<-no.tbl[,1];
    new.rows<-new.rows[!(rownames(new.rows) %in% rownames(meta)), , drop=FALSE];
    if (nrow(new.rows) > 0) {
      df<-data.frame(new.rows, stringsAsFactors=FALSE);
      colnames(df)<-colnames(meta);
      meta<-rbind(meta, df);
    }
  } 
  ####################################################################################
  
  saveRDS(meta, file=fn.meta);
  saveRDS(meta.lst, file=fn.meta.lst);
  
  meta;
}

# Retrieve a specified test statistics: p value, effect size, or allele frequency from downloaded dbGaP results files
RetrieveDbGapStat<-function(id.ana, id.std, stat.name=c('p value', 'effect size', 'allele frequency'), 
                            path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'),
                            own.table.min=20) {
  # id.ana, id.std        Matched analysis and study IDs to be retrieved
  # stat.name             The type of test statistics to summarize. Integer (1, 2, 3) or name of the statistics
  # path                  Path to input/output files  
  # own.table.min         Minimum number of analyses to make a study its own table
  
  own.table.min<-max(3, own.table.min); # At lease 3 analyses to make their own table
  
  # specify the type and name of test statistics to summarize
  stat.types<-c('p value' = 'phred', 'effect size' = 'es', 'allele frequency' = 'maf');
  stat.type<-stat.types[stat.name[1]];
  if (is.na(stat.type)) stat.type<-c('p value'=1);
  
  # Column names that fit to the type of test statistics
  cnm<-list(
    'phred' = c('p-value', 'p_value', 'p value', 'pvalue', 'p'),
    'es' = c('effect size', 'effect_size', 'effect-size', 'or', 'odds ratio', 'odds-ratio', 'odds_ratio', 'beta', 'gee', 'gee beta', 'gee-beta', 'gee_beta'),
    'maf' = c('maf', 'af', 'minor allele frequency', 'minor_allele_frequency', 'minor-allele-frequency', 'allele frequency', 'allele_frequency', 'allele-frequency')
  )[[stat.type]];
  
  
}


