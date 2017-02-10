##### Parse DisGeNET collection of disease-gene assocation

ParseDisGeNET<-function(url=c('Curated'="http://www.disgenet.org/ds/DisGeNET/results/curated_gene_disease_associations.tsv.gz", 
                              'BeFree'="http://www.disgenet.org/ds/DisGeNET/results/befree_gene_disease_associations.tsv.gz",
                              'All'="http://www.disgenet.org/ds/DisGeNET/results/all_gene_disease_associations.tsv.gz",
                              'PubMed'="http://www.disgenet.org/ds/DisGeNET/results/befree_results_only_version_4.0.tar.gz"),
                        path=paste(Sys.getenv("RCHIVE_HOME"), 'data/disease/public/disgenet', sep='/')) {
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  # Download source file
  fn.src <- sapply(strsplit(url, '/'), function(x) paste(path, 'src', x[length(x)], sep='/'));
  sapply(names(fn.src), function(nm) download.file(url[nm], fn.src[nm]))->x;
  
  # The full table
  tbls<-lapply(names(fn.src)[1:3], function(nm) {
    cat(nm, '\n');
    lns <- scan(fn.src[nm], what='', flush=TRUE, sep='\n', comment.char = '#');
    lns <- lns[grep('\t', lns)]; 
    cnm<-strsplit(lns[1], '\t')[[1]];
    tbl<-sapply(strsplit(lns[-1], '\t'), function(x) x[1:length(cnm)]);
    tbl[is.na(tbl)]<-'';
    tbl<-data.frame(t(tbl), stringsAsFactors = FALSE);
    names(tbl)<-cnm;
    cat("Loading in", nrow(tbl), 'rows\n');
    if ('score' %in% colnames(tbl)) tbl$score <- as.numeric(tbl$score); 
    if ('NofPmids' %in% colnames(tbl)) tbl$NofPmids <- as.numeric(tbl$NofPmids); 
    if ('NofSnps' %in% colnames(tbl)) tbl$NofSnps <- as.numeric(tbl$NofSnps); 
    tbl;
  }); 
  names(tbls) <- names(url)[1:3];
  fn <- sapply(names(tbls), function(nm) saveRDS(tbls[[nm]], paste(path, '/r/', tolower(nm), '.rds', sep='')))
  
  lns <- scan(fn.src[4], what='', flush=TRUE, sep='\n', comment.char = '#');
  lns <- lns[grep('\t', lns)]; 
  cnm <- strsplit(lns[1], '\t')[[1]];
  tbl <- sapply(strsplit(lns[-1], '\t'), function(x) x[1:length(cnm)]);
  tbl[is.na(tbl)] <- '';
  tbl <- data.frame(t(tbl), stringsAsFactors = FALSE);
  names(tbl) <- cnm; 
  if ('SECTION_NUM' %in% colnames(tbl)) tbl$SECTION_NUM <- as.numeric(tbl$SECTION_NUM); 
  if ('SENTENCE_NUM' %in% colnames(tbl)) tbl$SENTENCE_NUM <- as.numeric(tbl$SENTENCE_NUM); 
  saveRDS(tbl, paste(path, '/r/', tolower(names(url)[4]), '.rds', sep=''));
  
  mps <- lapply(tbls, function(tbl) split(tbl[, 'geneId'], tbl[, 'diseaseId'])); 
  mps[[4]] <- split(tbl$GENE_ID, tbl$DISEASE_ID); 
  names(mps[[4]]) <- paste('umls', names(mps[[4]]), sep=':'); 
  mps <- lapply(mps, function(m) lapply(m, unique)); 
  names(mps) <- names(url); 

  di <- do.call('rbind', lapply(tbls, function(t) t[, c('diseaseId', 'diseaseName')]));
  di <- di[!duplicated(di[, 1]), ];
  ns <- sapply(mps, function(m) sapply(m[di[, 1]], length)); 
  colnames(ns) <- paste('N', colnames(ns), sep='_'); 
  ns <- ns[, order(colMeans(ns))]; 
  di <- data.frame(row.names=di[, 1], stringsAsFactors = FALSE, Concept_ID=sub('umls:', '', di[, 1]), Name=di[, 2], ns);
  di$URL <- paste("http://www.ncbi.nlm.nih.gov/medgen", di[, 1], sep='/');
  di <- di[order(di[, 1]), ];

  mps <- mps[order(sapply(mps, function(x) length(unlist(x, use.names=FALSE))))]; 
  saveRDS(mps, file=paste(path, '/r/disease2gene_full', '.rds', sep=''));
  sapply(names(mps), function(nm) saveRDS(mps[[nm]], file=paste(path, '/r/disease2gene_', tolower(nm), '.rds', sep=''))) -> x; 
  saveRDS(di, file=paste(path, 'r', 'disease_summary.rds', sep='/'));
  
  di;
}

