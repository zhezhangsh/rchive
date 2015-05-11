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
      hd<-ln[substr(ln, 1, 2)=='# '];
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
  # id.ana, id.std        Matched analysis and study IDs, preferentially including full list of dbGaP analysis, both previously loaded and unloaded analysis
  # stat.name             The type of test statistics to summarize. Integer (1, 2, 3) or name of the statistics
  # path                  Path to input/output files  
  # own.table.min         Minimum number of analyses to make a study its own table
  
  library(rchive);
  
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
  )[stat.type];
  
  # Group analysis by study IDs
  if (length(id.ana) != length(id.std)) stop("Number of analyses and number of studies don't match!");
  names(id.ana)<-id.std;
  id.ana<-id.ana[!is.na(id.ana) & !is.na(id.std) & id.ana!='' & id.std!=''];
  std2ana<-split(as.vector(id.ana), id.std);
  std2ana<-lapply(std2ana, function(ids) {
    fn<-paste(path, '/r/full_table/', ids, '.rds', sep='');
    names(fn)<-ids;
    fn[file.exists(fn)];
  });
  std2ana<-std2ana[sapply(std2ana, length)>0];
  
  if (length(std2ana)>0) {
    # the table save the analyses of studies don't have their own table
    fn.tbl0<-paste(path, '/r/', names(cnm)[1], '_', 'phs000000.rds', sep='');
    if (file.exists(fn.tbl0)) tbl0<-readRDS(fn.tbl0) else tbl0<-matrix(nr=0, nc=0);
    ana0<-colnames(tbl0);
    
    fn.std<-paste(path, '/r/', names(cnm)[1], '_', names(std2ana), '.rds', sep=''); # file name of the tables storing stats of a study
    names(fn.std)<-names(std2ana);
    existing<-file.exists(fn.std); # whether the table corresponding to the study already exists
    std2ana.existing<-std2ana[existing]; # studies already have a table 
    std2ana.new<-std2ana[!existing]; # studies have no existing table 

    # add new analysis to existing tables
    if (length(std2ana.existing) > 0) {
      sapply(names(std2ana.existing), function(nm) {
        cat('Updating study', nm, '\n');
        tbl<-readRDS(fn.std[nm]);
        ana<-std2ana.existing[[nm]];
        ana<-ana[!(names(ana) %in% colnames(tbl))];
        if (length(ana) > 0) {
          tbl1<-LoadColumn(ana, cnm[[1]], row=c('snp id', 'snp_id', 'snp-id'), collapse='min', is.numeric=TRUE);
          tbl<-MergeTable(list(tbl, tbl1));
          saveRDS(tbl[!is.na(rowMeans(tbl, na.rm=TRUE)), !is.na(colMeans(tbl, na.rm=TRUE)), drop=FALSE], file=fn.std[nm]);
        }
      });
    }
    
    # new analysis to new tables
    if (length(std2ana.new) > 0) {
      tbl<-lapply(names(std2ana.new), function(nm) {
        cat('Updating study', nm, '\n');
        ana<-std2ana.new[[nm]];
        tbl.loaded<-tbl0[, colnames(tbl0) %in% names(ana), drop=FALSE];
        ana.unloaded<-ana[!(names(ana) %in% colnames(tbl0))];
        if (length(ana.unloaded) > 0) {
          tbl<-LoadColumn(std2ana.new[[nm]], cnm[[1]], row=c('snp id', 'snp_id', 'snp-id'), collapse='min', is.numeric=TRUE);
          tbl<-round(-10*log10(tbl));
          tbl<-MergeTable(list(tbl, tbl.loaded));
        } else tbl<-tbl.loaded;
        tbl;
      });
      
      # Save to new file or merge to table for orphan analyses
      for (i in 1:length(tbl)) {
        if (ncol(tbl[[i]]) >= own.table.min) {
          t<-tbl[[i]]; 
          t<-t[!is.na(rowMeans(t, na.rm=TRUE)), !is.na(colMeans(t, na.rm=TRUE)), drop=FALSE];
          saveRDS(t, fn.std[names(std2ana.new)[i]]);
        } else tbl0<-MergeTable(list(tbl0, tbl[[i]]));
      }
      
      # Remove analyses having their own table from orphan table
      not.orphan<-unlist(lapply(tbl, function(t) if (ncol(t) > own.table.min) colnames(t) else c()), use.names=FALSE);
      not.orphan<-unique(c(not.orphan, unlist(lapply(std2ana.existing, names), use.names=FALSE)));
      tbl0<-tbl0[, !(colnames(tbl0) %in% not.orphan), drop=FALSE];
      
      if (!setequal(ana0, colnames(tbl0))) saveRDS(tbl0[!is.na(rowMeans(tbl0, na.rm=TRUE)), !is.na(colMeans(tbl0, na.rm=TRUE)), drop=FALSE], file=fn.tbl0);
    }
  }
  
  std2ana;
}


