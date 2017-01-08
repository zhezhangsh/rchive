##### Parse DisGeNET collection of disease-gene assocation

ParseDisGeNET<-function(url=c('curated'="http://www.disgenet.org/ds/DisGeNET/results/curated_gene_disease_associations.txt.gz", 
                              'literature'="http://www.disgenet.org/ds/DisGeNET/results/literature_gene_disease_associations.txt.gz",
                              'befree'="http://www.disgenet.org/ds/DisGeNET/results/befree_gene_disease_associations.txt.gz",
                              'all'="http://www.disgenet.org/ds/DisGeNET/results/all_gene_disease_associations.txt.gz"), 
                        path=paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/disgenet', sep='/')) {
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  # Download source file
  fn.src<-sapply(strsplit(url, '/'), function(x) paste(path, 'src', x[length(x)], sep='/'));
  sapply(names(fn.src), function(nm) download.file(url[nm], fn.src[nm]))->x;
  
  # The full table
  tbls<-lapply(names(fn.src), function(nm) {
    cat(nm, '\n');
    lns<-scan(fn.src[nm], what='', flush=TRUE, sep='\n');
    cnm<-strsplit(lns[1], '\t')[[1]];
    tbl<-sapply(strsplit(lns[-1], '\t'), function(x) x[1:length(cnm)]); 
    tbl[is.na(tbl)]<-'';
    tbl<-data.frame(t(tbl), stringsAsFactors = FALSE);
    names(tbl)<-cnm;
    cat("Loading in", nrow(tbl), 'lines\n');
    tbl;
  }); 
  names(tbls)<-names(url);
  for (i in 1:length(tbls)) tbls[[i]][, 'score']<-as.numeric(tbls[[i]][, 'score']);
  fn<-sapply(names(tbls), function(nm) saveRDS(tbls[[nm]], paste(path, '/r/', nm, '.rds', sep='')))
  
  mp0<-split(tbls[[1]][, 'geneId'], tbls[[1]][, 'diseaseId']);
  mp1<-split(tbls[[4]][, 'geneId'], tbls[[4]][, 'diseaseId']);
  mp0<-lapply(mp0, unique);
  mp1<-lapply(mp1, unique);
  
  di<-do.call('rbind', lapply(tbls, function(t) t[, c('diseaseId', 'diseaseName')]));
  di<-di[!duplicated(di[, 1]), ];
  di<-data.frame(row.names=di[, 1], stringsAsFactors = FALSE, Concept_ID=sub('umls:', '', di[, 1]), Name=di[, 2],  
                 N_Curated=sapply(mp0[di[,1]], length), N_All=sapply(mp1[di[,1]], length));
  di$URL<-paste("http://www.ncbi.nlm.nih.gov/medgen", di[, 1], sep='/');
  di<-di[order(di[, 1]), ];
  
  tp<-strsplit(tbls[[4]][, 'associationType'], ', ');
  tp<-sort(unique(unlist(tp)));
  sr<-strsplit(tbls[[4]][, 'source'], ', ');
  sr<-sort(unique(unlist(sr)));
  
  t<-tbls[[1]];
  tp0<-sapply(tp, function(tp) grepl(tp, t[, 'associationType']));
  sr0<-sapply(sr, function(sr) grepl(sr, t[, 'source']));
  map0<-list(All=mp0);
  map0[(length(map0)+1):(length(map0)+length(tp))]<-lapply(tp, function(tp) split(t[tp0[, tp], 'geneId'], t[tp0[, tp], 'diseaseId']));
  map0[(length(map0)+1):(length(map0)+length(sr))]<-lapply(sr, function(sr) split(t[sr0[, sr], 'geneId'], t[sr0[, sr], 'diseaseId']));  
  names(map0)<-c('All', tp, sr);
  
  t<-tbls[[4]];
  tp1<-sapply(tp, function(tp) grepl(tp, t[, 'associationType']));
  sr1<-sapply(sr, function(sr) grepl(sr, t[, 'source']));
  map1<-list(All=mp1);
  map1[(length(map1)+1):(length(map1)+length(tp))]<-lapply(tp, function(tp) split(t[tp1[, tp], 'geneId'], t[tp1[, tp], 'diseaseId']));
  map1[(length(map1)+1):(length(map1)+length(sr))]<-lapply(sr, function(sr) split(t[sr1[, sr], 'geneId'], t[sr1[, sr], 'diseaseId']));  
  names(map1)<-c('All', tp, sr);
  
  saveRDS(map0, file=paste(path, 'r', 'disease2gene_curated.rds', sep='/'));
  saveRDS(map1, file=paste(path, 'r', 'disease2gene_all.rds', sep='/'));
  saveRDS(di, file=paste(path, 'r', 'disease_summary.rds', sep='/'));
  
  di;
}

