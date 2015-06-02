### Functions to query the SQLite database of dbSNP

# Query the dbSNP SQLite database by SNP IDs to get SNP position
GetSnpPosById<-function(ids, genome='GRCh38', path.db=paste(Sys.getenv("RCHIVE_HOME"), '/data/variant/db/variant_position.sqlite', sep='')) {

  library(dplyr);
  library(IRanges);
  library(GenomicRanges)
  
  db<-src_sqlite(path.db);
  tbl.nm<-src_tbls(db);
  
  genome<-tolower(genome);
  
  if (genome %in% c('hg20', 'hg38', 'grch38', 'human')) {
    t<-tbl(db, 'dbsnp_GRCh38');
  } else if (genome %in% c('hg19', 'grch37')) {
    t<-tbl(db, 'dbsnp_GRCh37');
  } else {
    warning('Unknown genome: ', genome, '\n');
    t<-c();
  }
  
  if (length(t) > 0) {
    t<-dplyr::filter(t, id %in% ids);
    t<-as.data.frame(collect(t));
    pos<-GRanges(t$chr, IRanges(start=t$pos, width=t$width));
    names(pos)<-as.vector(t$id);
  } else pos<-GRanges();
  
  pos;
}