#############################################################
# Parse NCBI Taxonomy
ParseTaxonomy<-function(ftp.file="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip", download.new=TRUE,
                        path=paste(Sys.getenv('RCHIVE_HOME'), 'data/taxonomy/public/ncbi', sep='/')) {

  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));
  
  # Download source file
  fn<-paste(path, 'src', 'ncbi_taxdmp.zip', sep='/');
  if (!file.exists(fn) | download.new) {
    download.file(ftp.file, fn);
    wd<-getwd();
    setwd(paste(path, 'src', sep='/'));
    try(unzip(fn));
    setwd(wd);
  }
  
  # read in species names
  nms.fn<-paste(path, 'src', 'names.dmp', sep='/');
  fll<-scan(nms.fn, sep='\n', what='', flush=TRUE);
  fld<-strsplit(fll, '\t');
  tbl0<-do.call(rbind, fld);
  tbl<-tbl0[, seq(1, 7, 2)];
  tbl[is.na(tbl)]<-'';
  tbl<-data.frame(tbl, stringsAsFactors=FALSE);
  names(tbl)<-c('ID', 'Name', 'Unique_Name', 'Class');
  
  # split
  id2name<-sapply(split(tbl[,2], tbl[,1]), unique);
  name2id<-sapply(split(tbl[,1], tbl[,2]), unique);
  cls<-split(tbl[, 1:3], tbl[[4]]);
  cmm<-cls[['common name']];
  id2name.cmm<-sapply(split(cmm[,2], cmm[, 1]), unique);
  name2id.cmm<-sapply(split(cmm[,1], cmm[, 2]), unique);
  
  # Representative names of each tax ID
  rep<-cls[["scientific name"]];
  id2name.rep<-rep[,2];
  names(id2name.rep)<-rep[,1];
  gbk<-sapply(split(cls[['genbank common name']][, 2], cls[['genbank common name']][, 1]), function(nm) nm[nchar(nm)==min(nchar(nm))][1]);
  cmm<-sapply(split(cls[['common name']][, 2], cls[['common name']][, 1]), function(nm) nm[nchar(nm)==min(nchar(nm))][1]);
  id2name.rep[names(cmm)]<-cmm;
  id2name.rep[names(gbk)]<-gbk;
  id2name.sht<-sapply(strsplit(id2name.rep, ' '), function(x) x[length(x)]);
  
  
  # Save data
  saveRDS(tbl, file=paste(path, 'r', 'full_table.rds', sep='/'));
  saveRDS(id2name, file=paste(path, 'r', 'id2name.rds', sep='/'));
  saveRDS(name2id, file=paste(path, 'r', 'name2id.rds', sep='/'));
  saveRDS(id2name.cmm, file=paste(path, 'r', 'id2name_common.rds', sep='/'));
  saveRDS(name2id.cmm, file=paste(path, 'r', 'name2id_common.rds', sep='/'));
  saveRDS(id2name.rep, file=paste(path, 'r', 'name_representative.rds', sep='/'));
  saveRDS(id2name.sht, file=paste(path, 'r', 'name_short.rds', sep='/'));
  
  saveRDS(cls, file=paste(path, 'r', 'table_by_class.rds', sep='/'));
  
  unique(tbl[,1]);
}

