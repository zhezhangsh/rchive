# Parse source files from dbSNP
#############################################################
ParseDbSNP<-function(nm, vcf, path=paste(RCHIVE_HOME, 'data/variant/public/dbsnp', sep='/')) {
  # nm              Name of a specific dbSNP database
  # vcf             Full path to the dbSNP vcf file on NCBI FTP server
  # path            home directory of data
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));
  
  library(VariantAnnotation);
  library(dplyr);
  library(DBI);
  
  fn<-paste(path, '/src/', nm, '_All.vcf.gz', sep='');
  fn0<-paste(path, '/r/full_vcf_', nm, '.rds', sep='');
  fn.pos<-paste(path, '/r/position_', nm, '.rds', sep='');

  # Download from source
  cat('Downloading ', nm, '...\n');
  download.file(vcf, fn);
  download.file(paste(vcf, '.tbi', sep=''), paste(fn, '.tbi', sep=''));
  
  # load VCF
  cat('Loading ', nm, '...\n'); 
  dbsnp<-scanVcf(fn, param=ScanVcfParam(fixed='ALT', info=NA, geno=NA))[[1]];
  saveRDS(dbsnp, file=fn0);

  # SNP positions
  gr<-dbsnp$rowData;
  library(GenomicRanges);
  if (file.exists(fn.pos))  {
    id.old<-names(readRDS(fn.pos));
    log0<-list(N=length(gr), Source=vcf, Added=setdiff(names(gr), id.old), Removed=setdiff(id.old, names(gr))); # Get update info.
  } else {
    log0<-list(N=length(gr), Source=vcf, Added=names(gr), Removed=c()); # Get update info.
  }
  saveRDS(gr, file=fn.pos);
  
  # save update log
  fn.log<-paste(path, 'log.rds', sep='/');
  if (file.exists(fn.log)) log<-readRDS(fn.log) else log<-list();
  log[[as.character(nm)]][[as.character(Sys.Date())]]<-log0;
  saveRDS(log, fn.log);
  
  # create database table
  # create SQLite database or connect to existing one
  db.fn<-paste(RCHIVE_HOME, '/data/variant/db/variant_position.sqlite', sep='');
  if (!file.exists(db.fn)) db<-src_sqlite(db.fn, create=TRUE) else db<-src_sqlite(db.fn);  
  cat('Creating table', nm, '...\n');
  tbl<-data.frame(id=names(gr), chr=as.vector(seqnames(gr)), pos=start(gr), width=width(gr), stringsAsFactors=FALSE);
  tbl.nm<-paste('dbsnp', nm, sep='_');
  tbls<-src_tbls(db);
  if (tbl.nm %in% tbls) dbRemoveTable(db$con, tbl.nm);
  t<-copy_to(db, tbl, temporary=FALSE, name=, indexes=list('id', 'chr', 'pos'));
}

#############################################################
