require(rchive);
require(RCurl); 

# Download and parse the MSigDb database
path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/regulatorycircuits', sep='/');
dir.create(path, showWarnings=FALSE);
dir.create(paste(path, 'r', sep='/'), showWarnings=FALSE);
dir.create(paste(path, 'src', sep='/'), showWarnings=FALSE);

