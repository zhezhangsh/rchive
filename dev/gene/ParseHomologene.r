# Parse the NCBI Homologene genes of multiple species
ParseHomologene<-function(ftp.file="ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",
                          path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/homologene', sep='/')) {
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  fn.ftp<-paste(path, 'src', 'homologene.data', sep='/');
  download.file(ftp.file, fn.ftp);
  
  homolo<-scan(fn.ftp, sep='\n', flush=TRUE, what='');
  homolo<-strsplit(homolo, '\t');
  homolo<-do.call('rbind', homolo);
  homolo<-data.frame(homolo, stringsAsFactors=FALSE);
  names(homolo)<-c('HID', 'Taxonomy_ID', 'Gene_ID', 'Gene_Symbol', 'Protein_GI', 'Protein_ACC');
  saveRDS(homolo, file=paste(path, 'r', 'homologene.rds', sep='/'));
  
  # Required files
  fn.tax<-paste(Sys.getenv("RCHIVE_HOME"), 'data/taxonomy/public/ncbi/r/name_short.rds', sep='/');
}