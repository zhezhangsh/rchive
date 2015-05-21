################################################################################
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
    meta<-data.frame(matrix('', nr=0, nc=7, dimnames=list(c(), c('name', 'study', 'genome', 'dbsnp', 'description', 'method', 'file'))), stringsAsFactors=FALSE);
  }
  # Original metadata with each ftp file of analysis; as a list
  fn.meta.lst<-paste(path, 'r', 'original_metadata_list.rds', sep='/');
  if (file.exists(fn.meta.lst)) meta.lst<-readRDS(fn.meta.lst) else {
    meta.lst<-list();
  }
  
  # Download and parse metadata table
  ln<-strsplit(getURL(url), '\n')[[1]];
  gw<-do.call('rbind', strsplit(ln, '\t'));
  colnames(gw)<-gw[1,];
  gw<-gw[-1, , drop=FALSE];
  gw<-data.frame(gw, stringsAsFactors=FALSE);
  gw<-gw[!(as.numeric(gw[,1]) %in% as.numeric(redundant.studies)), ];
  
  # Each analysis in the current meta table
  for (i in 1:nrow(gw)) {
    id<-rownames(gw)[i];
    fn.ftp<-gw[id, 'analysis.result']; # Full path to FTP URL
    fn.gz<-strsplit(fn.ftp, '/')[[1]]; # file name
    id.study<-strsplit(fn.gz[length(fn.gz)], '\\.')[[1]][1];
    id.analysis<-strsplit(fn.gz[length(fn.gz)], '\\.')[[1]][2];    
    fn.gz<-paste(path, 'src', fn.gz[length(fn.gz)], sep='/'); # local copy of the analysis source file
    
    fn.r<-paste(path, '/r/full_table/', id.analysis, '.rds', sep=''); # analysis file loaded into R
    
    # Whether the analysis has been fully loaded
    loaded<-file.exists(fn.gz) & file.exists(fn.r) & id.analysis %in% rownames(meta) & id.analysis %in% names(meta.lst);
    
    # Load the analysis
    if (!loaded | update.all.analysis) {
      cat('Loading analysis', id.analysis, '\n');
      if (!file.exists(fn.gz) | update.all.analysis) download.file(fn.ftp, fn.gz);
      
      # load in table and save as R
      ln<-scan(fn.gz, what='', flush=TRUE, sep='\n');
      # very tediously parse file header and add header info to metadata table and list
      hd<-ln[substr(ln, 1, 2)=='# '];
      hd<-sub('^# ', '', hd);
      hd<-strsplit(hd, ':\t');
      hd<-lapply(hd, function(hd) sub('^ ', '', hd));
      hd<-lapply(hd, function(hd) sub(' $', '', hd));      
      hd<-do.call('rbind', hd);
      hd.nm<-tolower(hd[,1]);
      hd<-hd[,2];
      names(hd)<-hd.nm;
      
      # Save table
      cnm<-strsplit(ln[length(hd)+1], '\t')[[1]];
      tbl<-strsplit(ln[(length(hd)+2):length(ln)], '\t');
      tbl<-t(sapply(tbl, function(x) x[1:length(cnm)]));
      tbl<-data.frame(tbl, stringsAsFactors=FALSE); # all columns are character
      colnames(tbl)<-cnm;
      saveRDS(tbl, file=fn.r);
      
      # Add original metadata
      meta.lst[[id.analysis]]<-hd;
      
      # Add new row or update a row
      row<-as.vector(c(hd['name'], id.study, hd[c('human genome build', 'dbsnp build', 'description', 'method')], fn.ftp));
      if (id.analysis %in% rownames(meta)) meta[id.analysis, ]<-row else {
        meta<-rbind(meta, row);
        rownames(meta)[nrow(meta)]<-id.analysis;
      }
    }
  }

  saveRDS(meta[order(rownames(meta)), ], file=fn.meta);
  saveRDS(meta.lst[order(names(meta.lst))], file=fn.meta.lst);
  
  meta;
}

