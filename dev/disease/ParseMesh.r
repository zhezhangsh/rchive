#### Parse MeSH diseases and chemicals from text files

############################################
# Parse MeSH diseases
ParseMeshDisease<-function(ftp="ftp://nlmpubs.nlm.nih.gov/online/mesh/.asciimesh/d2015.bin",
                           path=paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/mesh', sep='/')) {
  mesh.d<-ParseMesh(ftp, path);
  mesh<-data.frame(row.names=names(mesh.d), 
                   Name = sapply(mesh.d, function(x) x['MH'][[1]]), 
                   Type = sapply(mesh.d, function(x) x['RECTYPE'][[1]]),
                   Date = sapply(mesh.d, function(x) x['DA'][[1]]),
                   URL = paste("http://www.ncbi.nlm.nih.gov/mesh/?term=", names(mesh.d), sep=''),
                   Description = sapply(mesh.d, function(x) {desc<-x[['MS']]; if (length(desc)==0) '' else desc[[1]]; }),
                   stringsAsFactors=FALSE);
  mesh[is.na(mesh)]<-'';
  saveRDS(mesh, file=paste(path, 'r', 'mesh.rds', sep='/'));
  mesh;
}

############################################
# Parse MeSH chemicals
ParseMeshChemical<-function(ftp="ftp://nlmpubs.nlm.nih.gov/online/mesh/.asciimesh/c2015.bin", 
                            path=paste(Sys.getenv("RCHIVE_HOME"), 'data/chemical/public/mesh', sep='/')) {
  mesh.c<-ParseMesh(ftp, path);
  mesh<-data.frame(row.names=names(mesh.c), 
                   Name = sapply(mesh.c, function(x) x['NM'][[1]]), 
                   Type = sapply(mesh.c, function(x) x['RECTYPE'][[1]]),
                   Date = sapply(mesh.c, function(x) x['DA'][[1]]),
                   URL = paste("http://www.ncbi.nlm.nih.gov/mesh/?term=", names(mesh.c), sep=''),
                   Description = sapply(mesh.c, function(x) {desc<-x[['NO']]; if (length(desc)==0) '' else desc[[1]]; }),
                   stringsAsFactors=FALSE);
  mesh[is.na(mesh)]<-'';
  saveRDS(mesh, file=paste(path, 'r', 'mesh.rds', sep='/'));
  mesh;
}

############################################
ParseMesh<-function(ftp, path) {
  # ftp     URL to source FTP file
  # path    Path to output files
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  # Download and load file
  fn.ftp<-paste(path, '/src/', sapply(strsplit(ftp, '/'), function(x) x[length(x)]), sep='');
  download.file(ftp, fn.ftp);
  lns<-scan(fn.ftp, what='', sep='\n', flush=TRUE);
  
  # Split by disease
  ind<-which(lns == '*NEWRECORD');
  ids<-rep(1:length(ind), c(ind[-1]-ind[-length(ind)]-1, length(lns)-ind[length(ind)]));
  fld<-strsplit(lns[-ind], ' = ');
  val<-sapply(fld, function(x) x[2]);
  names(val)<-sapply(fld, function(x) x[1]);
  names(val)[names(val)=='UI']<-'ID';
  grp<-split(val, ids);
  grp<-lapply(grp, function(grp) split(as.vector(grp), names(grp)));
  names(grp)<-sapply(grp, function(grp) grp['ID'][[1]]);
  
  saveRDS(grp, file=paste(path, 'r', 'mesh_by_id.rds', sep='/'));
  
  grp;
}