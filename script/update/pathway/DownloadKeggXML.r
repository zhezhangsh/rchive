path.out<-paste(Sys.getenv('RCHIVE_HOME'), 'data/pathway/public/kegg/src', sep='/');
if (!file.exists(path.out)) dir.create(path.out, recursive = TRUE); 

library(gage);
library(pathview);

data(bods);
sps<-bods[, 'kegg code'];
ids<-lapply(sps, function(s) names(gage::kegg.gsets(s)[[1]]));
ids<-lapply(ids, function(ids) sapply(strsplit(ids, ' '), function(x) x[1]));
names(ids)<-sps;

status<-lapply(sps, function(s) {download.kegg(sub(s, '', ids[[s]]), s, path.out)})

