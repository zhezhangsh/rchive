################################################################################
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
}
