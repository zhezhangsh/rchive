# Parse the NCBI Unigene genes of one species
ParseUnigene<-function(ftp.file, species='human', download.new=FALSE, 
                       col.names=c('ID', 'GENE_ID', 'GENE', 'CHROMOSOME', 'CYTOBAND', 'TITLE'),
                       path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/unigene', sep='/')) {

  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  # Download FTP file
  fn<-paste(path, 'src', paste(species, '.data.gz', sep=''), sep='/');   
  if (download.new | !file.exists(fn)) download.file(ftp.file, fn); 
  
  cat("Loading source file...\n");
  lns<-scan(fn, what='', sep='\n', flush=TRUE);   
  ind<-grep('//', lns); # break lines between genes
  nln<-ind-c(0, ind[-length(ind)]); # Number of lines per gene
  ign<-rep(1:length(nln), nln); # Gene index of each line
  
  # Split into name/value pairs
  spc<-regexpr(' ', lns); 
  spc[ind]<-2; 
  fld<-substr(lns, 1, spc-1); # Field name
  val<-substr(lns, spc+1, nchar(lns)); # Field value
  val<-stringr::str_trim(val, 'left'); # remove whitespaces at the beginning
  names(fld)<-names(val)<-ign;

  #### Make the gene annotation table
  cat('Creating gene annotation table ...\n'); 
  if (!('ID' %in% col.names)) col.names<-c('ID', col.names); 
  cls<-lapply(col.names, function(cnm) {
    v<-val[fld==cnm]; 
    split(as.vector(v), names(v)); 
  }); 
  names(cls)<-col.names;
  ids<-cls[['ID']];
  ids<-sapply(ids, function(id) id[1]); 
  ids<-ids[order(as.numeric(names(ids)))];
  tmp<-rep('', length(ids));
  names(tmp)<-names(ids); 
  tbl<-sapply(cls[col.names!='ID'], function(cls) {
    c<-paste(cls, sep=';'); 
    names(c)<-names(cls); 
    c<-c[names(c) %in% names(tmp)];
    tmp[names(c)]<-c; 
    as.vector(tmp); 
  }); 
  rownames(tbl)<-ids;
  tbl<-data.frame(tbl, stringsAsFactors = FALSE); 
  saveRDS(tbl, file=paste(path, 'r', paste(species, '_genes.rds', sep=''), sep='/')); 
  out<-list(annotation=tbl); 
  
  # Tissue specific expression
  cat('Create tissue expression table ...\n');
  expr<-val[fld=='EXPRESS'];
  if (length(expr) > 0) {
    expr<-strsplit(expr, '\\|'); 
    expr.ind<-rep(names(expr), sapply(expr, length)); 
    expr.val<-unlist(expr, use.names=FALSE);
    expr.val<-stringr::str_trim(expr.val); 
    expr.cnm<-sort(unique(expr.val));
    expr.tbl<-matrix(0, nr=length(ids), nc=length(expr.cnm), dimnames=list(names(ids), expr.cnm)); 
    for (i in 1:length(expr.cnm)) {
      j<-unique(expr.ind[expr.val==expr.cnm[i]]); 
      expr.tbl[j, i]<-1;
    };
    rownames(expr.tbl)<-ids;
    saveRDS(expr.tbl, file=paste(path, 'r', paste(species, '_expression.rds', sep=''), sep='/')); 
    out$expression<-expr.tbl;
  }
 
  # Protein similarity to other species
  cat('Creating gene homolog table ...\n'); 
  prt<-val[fld=='PROTSIM']; 
  if (length(prt) > 0) {
    prt<-strsplit(prt, '; '); 
    prt<-do.call('rbind', prt); 
    prt.cnm<-c('ORG', 'PROTGI', 'PROTID', 'PCT', 'ALN'); 
    for (i in 1:ncol(prt)) prt[, i]<-sub(paste(prt.cnm[i], '=', sep=''), '', prt[, i]); 
    tbl.prt<-cbind(ids[rownames(prt)], prt); 
    rownames(tbl.prt)<-1:nrow(tbl.prt);
    tbl.prt<-data.frame(tbl.prt, stringsAsFactors = FALSE); 
    names(tbl.prt)<-c('ID', prt.cnm);
    tbl.prt$PCT<-as.numeric(tbl.prt$PCT); 
    tbl.prt$ALN<-as.numeric(tbl.prt$ALN); 
    saveRDS(tbl.prt, file=paste(path, 'r', paste(species, '_homolog.rds', sep=''), sep='/')); 
    out$homolog=tbl.prt;
  }
  
  out; 
}