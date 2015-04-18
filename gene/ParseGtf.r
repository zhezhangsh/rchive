# Load a GTF file to create a GRanges
ParseGtf<-function(fn.in, fn.out, path=paste(RCHIVE_HOME, 'data/gene/public/gtf/r', sep='/'), Dbxref=TRUE, group=TRUE, separators=c()) {
  # fn.in, fn.out		Name of input/output files
  # path    			  Path to output file
  # Dbxref, group		Whether to further parse these attribute fields
  # separators  		The separators to be used to parse Dbxref or group field
  
  library(GenomicRanges);
  library(rtracklayer);
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  
  cat('Importing GTF file:', fn.in, '\n');
  gr<-rtracklayer::import(fn.in);
  cat('The loaded GTF file includes', length(gr), 'entries.\n');
  names(gr)<-1:length(gr);
  
  meta<-elementMetadata(gr);
  rownames(meta)<-names(gr);
  
  ##############################################################################
  # Function that performs the splitting
  splitColumn<-function(splt, separator) {
    n<-sapply(splt, length);
    pairs<-strsplit(unlist(splt, use.names=FALSE), separator);
    nm<-sapply(pairs, function(x) x[1]);
    cnm<-unique(nm);
    vl<-sapply(pairs, function(x) x[length(x)]);
    names(vl)<-rep(names(splt), n);
    mp<-split(vl, nm);
    df<-matrix('', nrow=length(gr), ncol=length(cnm), dimnames=list(names(gr), cnm));
    for (i in 1:length(cnm)) df[names(mp[[cnm[i]]]), cnm[i]]<-as.vector(mp[[cnm[i]]]); 
    df<-data.frame(df, stringsAsFactors=FALSE);
    colnames(df)<-tolower(colnames(df));
    df;
  }
  ##############################################################################
  
  # Further parse 'Dbxref' column
  if (Dbxref) {
    dbx<-as.vector(meta$Dbxref);
    if (length(separators) !=2) separators<-c('', ':');
    names(dbx)<-names(gr);
    df<-splitColumn(dbx, separators[2]);
    meta<-cbind(meta[, !(colnames(meta) %in% 'Dbxref')], df);
  }
  
  # Further parse 'group' column	
  if (group & 'group' %in% colnames(meta)) {
    grp<-as.vector(meta$group);
    if (length(separators) !=2) separators<-c('; ', ' ');
    splt<-strsplit(grp, separators[1]);
    names(splt)<-names(gr);
    df<-splitColumn(splt, separators[2]);
    meta<-cbind(meta[, !(colnames(meta) %in% 'group')], df);
  }
  
  elementMetadata(gr)<-meta;
  saveRDS(gr, file=paste(path, '/', fn.out, '_full.rds', sep=''));
  
  # Minimal metadata for slim version
  cnm<-strsplit("source;type;score;phase;transcript_id;transcript_name;transcript_type;gene_id;gene_name;gene_type", ';')[[1]];
  colnames(meta)<-tolower(colnames(meta));
  meta<-meta[, colnames(meta) %in% cnm];
  cnm<-cnm[cnm %in% colnames(meta)];
  meta<-meta[, cnm];
  elementMetadata(gr)<-meta;
  saveRDS(gr, file=paste(path, '/', fn.out, '_slim.rds', sep=''));
  
  # split by type
  if ('type' %in% colnames(meta)) {
    tp<-as.vector(meta[['type']]);
    meta<-meta[, colnames(meta)!='type'];
    elementMetadata(gr)<-meta;
    gr0<-split(gr, tp);
    f<-sapply(names(gr0), function(nm) saveRDS(gr0[[nm]], file=paste(path, '/', fn.out, '_', tolower(nm), '.rds', sep='')));
  }
  
}