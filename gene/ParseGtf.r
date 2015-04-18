# Load a GTF file to create a GRanges
ParseGtf<-function(fn.in, fn.out, path.out='.', types=c('exon', 'gene', 'transcript'), Dbxref=TRUE, group=TRUE, separators=c()) {
  # fn.in, fn.out		Name of input/output files
  # path.out			  Path to output file
  # types				    The types of regions to be saved separately
  # Dbxref, group		Whether to further parse these attribute fields
  # separators		The separators to be used to parse Dbxref or group field
  
  library(GenomicRanges);
  library(rtracklayer);
  
  if (!file.exists(path.out)) dir.create(path.out, recursive=TRUE);
  
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
  
  saveRDS(gr, file=paste(path.out, '/', fn.out, '_full.rds', sep=''));
  if (length(type)>0 & grepl('type', colnames(meta))) {
    grs<-split(gr, as.vector(gr$type));
    cat('Saving', length(grs), 'region types to output folder.\n');
    fn<-sapply(names(grs), function(nm) saveRDS(grs[[nm]], file=paste(path.out, '/', fn.out, '_', nm, '.rdata', sep='')));
  }
}