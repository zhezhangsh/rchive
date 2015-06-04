# Parse CHEBI chemical records

ParseChebi<-function(ftp="ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology//chebi.obo", 
                     path=paste(Sys.getenv("RCHIVE_HOME"), 'data/chemical/public/chebi', sep='/')) {
  # ftp     URL to source FTP file
  # path    Path to output files
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  # Download the .obo ontology file
  fn<-sapply(strsplit(ftp, '/'), function(x) x[length(x)]); 
  fn.loc<-paste(path, 'src', fn, sep='/');
  download.file(ftp, fn.loc);
  
  # load and parse the .obo data
  lns<-scan(fn.loc, what='', sep='\n', flush=TRUE); 
  ind<-which(lns == '[Term]');
  
  ent<-rep(0:length(ind), c(ind[1]-1, ind[-1]-ind[-length(ind)], length(lns)+1-ind[length(ind)])); # entry ID
  ent[ind]<-0;
  ent.id<-as.character(1:length(ind));
  
  fld<-strsplit(lns, ": "); # separate key:value pair
  key<-sapply(fld, function(x) x[1]);
  val<-sapply(fld, function(x) x[2]);
  names(val)<-ent;
  # split by key types
  by.typ<-split(val[ent!=0], key[ent!=0]);
  
  # make the ChEBI table
  id<-sapply(by.typ[['id']], function(x) x[1])[ent.id];
  name<-sapply(by.typ[['name']], function(x) x[1])[ent.id];
  def<-sapply(by.typ[['def']], function(x) sapply(strsplit(x, '\"'), function(x) x[2]))[ent.id];
  def[is.na(def)]<-'';
  url<-paste("http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", id);
  chebi<-data.frame(row.names=id, Name=name, URL=url, Definition=def, stringsAsFactors=FALSE);
  saveRDS(chebi, file=paste(path, 'r', 'chebi.rds', sep='/'));
  
  # synonym table
  syn1<-sub('^\"', '', sapply(strsplit(by.typ[['synonym']], '\\" '), function(x) x[1]));
  syn2<-sapply(strsplit(by.typ[['synonym']], '\\" '), function(x) x[2]);
  syn3<-sapply(strsplit(syn2, ' \\['), function(x) x[1]);
  syn4<-sub(':\\]', '', sapply(strsplit(syn2, ' \\['), function(x) x[2]));
  syn<-data.frame(ID=id[names(syn1)], Synonym=syn1, Relationship=syn3, Source=syn4, stringsAsFactors=FALSE);
  saveRDS(syn, file=paste(path, 'r', 'chebi_synonym.rds', sep='/'));
  
  # alternative ID
  alt<-by.typ[['alt_id']];
  alt2id<-id[names(alt)];
  names(alt2id)<-alt;
  saveRDS(alt2id, file=paste(path, 'r', 'chebi_alternative_id.rds', sep='/'));
  
  # list by ID
  names(val)<-key;
  by.id<-split(val[ent>0], ent[ent>0]);
  names(by.id)<-id[names(by.id)];
  saveRDS(by.id, file=paste(path, 'r', 'chebi_by_id.rds', sep='/'));
  
  chebi;
}