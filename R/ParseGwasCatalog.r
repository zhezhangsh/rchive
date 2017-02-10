# Download and process GWAS results from GWAS catalog
ParseGwasCatalog<-function(url="https://www.ebi.ac.uk/gwas/api/search/downloads/full", 
                           path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/gwascatalog', sep='/')) {
  # url    Source file URL
  # path   Path to output files
  # update.all.pubmed   Re-download all PubMed entries if TRUE
  
  options()
  
  library(RCurl);
  library(rchive);
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  # download the table 
  ln <- strsplit(getURL(url), '\n')[[1]];
  gw <- do.call('rbind', strsplit(ln, '\t'));
  colnames(gw)<-gw[1,];
  gw <- gw[-1, ];
  gw <- data.frame(gw, stringsAsFactors=FALSE);
  
  # this will download an outdated version of GWAScatalog, 
  # but its column names were used across the rest of this function
  # so it's downloaded just to provide the column names
  # library(NCBI2R);
  # gw0 <- GetPublishedGWAS(); 
  # if (ncol(gw) == ncol(gw0)) colnames(gw) <- colnames(gw0) else warning("Unmatched column names between versions!")
  
  # save a copy of full table
  fn.src <- sapply(strsplit(url, '://'), function(x) x[2]);
  fn.src <- gsub('/', '_', fn.src);
  
  ##########################################################################################################
  # start working on the table
  
  # fields of numeric values
  numeric.fields <- c('CHR_POS', 'MERGED', 'INTERGENIC', 'P.VALUE', 'PVALUE_MLOG', 
                      'OR.or.BETA', 'UPSTREAM_GENE_DISTANCE', 'DOWNSTREAM_GENE_DISTANCE'); # fields shold be numeric
  numeric.fields <- numeric.fields[numeric.fields %in% colnames(gw)]; 
  num <- gw[, numeric.fields];
  num <- apply(num, 2, as.numeric);
  gw[, numeric.fields] <- num;
  colnames(gw) <- gsub('\\.', '_', colnames(gw));
  colnames(gw) <- gsub('__', '_', colnames(gw));
  colnames(gw) <- gsub('^_', '', colnames(gw));
  colnames(gw) <- gsub('_$', '', colnames(gw));
  colnames(gw) <- gsub('^X', '', colnames(gw));
  
  # fields must have valid value
  required.fields<-c('SNP_ID_CURRENT', 'CHR_ID', 'CHR_POS', 'PVALUE_MLOG');
  req<-gw[, required.fields];
  n<-apply(req, 2, function(x) is.na(x) | x=='');
  # cnv<-gw[gw$CNV=='Y', , drop=FALSE]; # CNV variants
  inc <- gw[rowSums(n)>0, , drop=FALSE];
  gw0 <- gw; 
  gw  <- gw0[gw$CNV!='Y' & rowSums(n)==0, , drop=FALSE]; 
  
  saveRDS(gw0, file=paste(fn.src, '.rds', sep=''));
  # saveRDS(cnv, file=paste(path, 'r', 'gwascatalog_just_cnv.rds', sep='/'));
  saveRDS(inc, file=paste(path, 'r', 'gwascatalog_incomplete_info.rds', sep='/'));
  saveRDS(gw, file=paste(path, 'r', 'gwascatalog_complete_info.rds', sep='/'));
  
  #####################
  gw<-data.frame(id=paste('rs', gw[['SNP_ID_CURRENT']], sep=''), phred=as.integer(round(10*abs(gw[['PVALUE_MLOG']]))), gw, stringsAsFactors=FALSE);
  #####################
  
  # the combination of SNP ID, PubMed ID, Trait, and p value text should be unique
  dup.fields<-c('REPORTED_GENE_S', 'MAPPED_GENE', 'REGION', 'SNP_GENE_IDS'); # the fields that might cause duplicated entries
  gw<-gw[!duplicated(gw[, !(colnames(gw) %in% dup.fields)]), ];
  uniq<-paste(gw[['id']], gw[['PUBMEDID']], gw[['DISEASE_TRAIT']], gw[['P_VALUE_TEXT']], sep='@_@');
  dup<-uniq[duplicated(uniq)];# duplicated entries
  dups<-lapply(dup, function(dup) gw[uniq==dup, ]);
  dups<-sapply(dups, function(d) rownames(d[order(d[, 'phred'])])[nrow(d)]); # pick one with the highest phred score
  gw<-rbind(gw[!(uniq %in% dup), ], gw[dups,]); # duplicated enries removed
  gw<-gw[!is.na(gw$phred), ]; # Entry will be removed if reported p value is 0
  
  ###############################################################################################################
  ###############################################################################################################
  cat('Re-formatting data ...\n');
  
  # re-format p value text
  txt<-gw[['P_VALUE_TEXT']];
  txt[is.na(txt)]<-'';
  txt<-gsub(' ', '', txt);
  txt<-gsub("[()]", '', txt);
  #anno[,4]<-txt;
  
  # unique analysis and study
  ana<-paste(gw[['PUBMEDID']], gw[['DISEASE_TRAIT']], gw[['P_VALUE_TEXT']], sep='@_@'); # analyses
  sty<-split(ana, gw[['PUBMEDID']]); # studies
  sty<-lapply(sty, unique);
  sty<-lapply(sty, sort);
  
  # analysis ID
  ana.nm<-unlist(sty, use.names=FALSE); # analysis names
  ana.id<-paste('pmid', rep(names(sty), sapply(sty, length)), '_', unlist(sapply(sty, function(x) 1:length(x))), sep=''); # analysis ID
  names(ana.id)<-ana.nm;
  id<-as.vector(ana.id[ana]); # analysis ID of each row
  
  # SNP positions
  pos<-gw[, c('id', 'CHR_ID', 'CHR_POS')];
  pos<-pos[!duplicated(pos[[1]]), ];
  pos<-pos[order(pos[,2], pos[,3]), ];
  
  # output table
  tbl<-matrix(nr=length(unique(gw$id)), nc=length(ana.id), dimnames=list(as.vector(pos[,1]), sort(ana.id)));
  
  # phred score of variants in each analysis
  phred<-gw$phred;
  names(phred)<-as.vector(gw$id);
  phred<-split(phred, id);
  for (i in 1:length(phred)) {
    v<-phred[[i]];
    p<-split(as.vector(v), names(v));
    p<-sapply(p, max); 
    tbl[names(p), names(phred)[i]]<-as.integer(p);
    #print(i); 
  }
  
  # Phred score matrix
  pos[[3]]<-as.integer(pos[[3]]);
  tbl<-data.frame(pos, tbl);
  rownames(tbl)<-as.vector(tbl[[1]]);
  #colnames(tbl)[1:3]<-c('id', 'chr', 'pos');
  tbl<-tbl[, 4:ncol(tbl)];
  saveRDS(tbl, file=paste(path, 'r/phred_gwascatalog.rds', sep='/'));
  
  ###############################################################################################################
  ###############################################################################################################
  cat('Getting study information from Pubmed ...\n');
  
  # annotate analyses
  # an<-gw[!duplicated(gw$PUBMEDID), c(4, 7, 6, 8, 35, 9, 10, 11, 12)]; # select fields
  an<-gw[!duplicated(id), ]; # one analysis, one row
  rownames(an)<-id[!duplicated(id)];
  sty.id<-paste('pmid', an$PUBMEDID, sep=''); # Non-unique study ID
  st.id<-unique(sty.id); # unique study ID
  id<-id0<-sub('^pmid', '', st.id);
  
  library(annotate); 
  if (file.exists(paste(path, 'r/pubmed_downloaded.rds', sep='/'))) {
    pm0 <- readRDS(paste(path, 'r/pubmed_downloaded.rds', sep='/'));
    id0 <- id[!(id %in% names(pm0))];
  } else {
    pm0 <- list();
  }; 
  pm <- list();
  if (length(id0)>0) {
    for (i in 1:length(id0)) {
      pub <- pubmed(id0[i], disp='data'); 
      pm[[i]] <- getPMInfo(pub)[[1]]; 
    };    
    names(pm) <- id0; 
    pm0 <- c(pm0, pm); 
    saveRDS(pm0, paste(path, 'r/pubmed_downloaded.rds', sep='/')); 
  }; 

  # PubMed table
  ttl <- sapply(pm0, function(x) x$title);
  jnl <- sapply(pm0, function(x) x$MedlineTA);
  dat <- sapply(pm0, function(x) x$JrnlInfo['year']);
  url <- paste("http://www.ncbi.nlm.nih.gov/pubmed/?term=", names(pm0), sep='');
  abs <- sapply(pm0, function(x) {a <- x$abstract; if (is.null(a)) a <- ''; a;});
  pubmed <- data.frame(Title=ttl, Journal=jnl, Date=dat, URL=url, Abstract=abs, stringsAsFactors = FALSE); 
  saveRDS(pubmed, file=paste(path, 'r/pubmed.rds', sep='/'));
  
  id2pub<-lapply(rownames(pubmed), function(id) c(ID=id, pubmed[id, ]));
  names(id2pub)<-rownames(pubmed);
  saveRDS(id2pub, file=paste(path, 'r/pubmed_by_id.rds', sep='/'));
  
  ###############################################################################################################
  ###############################################################################################################
  cat('Putting together analysis and study tables ...\n');
  
  # Number of phred scores of each analysis
  n<-sapply(rownames(an), function(id) {
    p<-tbl[, id];
    length(p[!is.na(p)]);
  });
  mx<-apply(tbl, 2, function(x) max(x, na.rm=TRUE));
  mn<-apply(tbl, 2, function(x) min(x, na.rm=TRUE));
  mx<-sapply(rownames(an), function(id) max(tbl[, id], na.rm=TRUE));
  mn<-sapply(rownames(an), function(id) min(tbl[, id], na.rm=TRUE));
  
  # analysis annotation table
  title<-pubmed[, 'Title'];
  names(title)<-paste('pmid',rownames(pubmed), sep='');
  desc<-paste(as.vector(an$Study), 'Disease trait:', as.vector(an$DiseaseTrait));
  ana<-data.frame(Source=rep('GWAScatalog', nrow(an)), Name=as.vector(title[sty.id]), Study_ID=sty.id,
                  URL=paste('http://www.ncbi.nlm.nih.gov/pubmed/', sub('^pmid', '', sty.id), sep=''),
                  Num_Variants=as.vector(n), Min_Phred=as.vector(mn), Max_Phred=as.vector(mx), Description=desc, row.names=rownames(an));
  ana<-ana[order(rownames(ana)), , drop=FALSE];
  saveRDS(ana, file=paste(path, 'r/analysis.rds', sep='/'))
  
  # study anotation table
  n.ana<-as.vector(table(ana$Study)[paste('pmid', rownames(pubmed), sep='')]);
  n.ana[is.na(n.ana)]<-0;
  std<-data.frame(Source=rep('GWAScatalog', nrow(pubmed)), Name=as.vector(pubmed$Title), Num_Analyses=n.ana, URL=as.vector(pubmed$URL), 
                  Description=as.vector(pubmed$Abstract), row.names=paste('pmid', rownames(pubmed), sep=''), stringsAsFactors=FALSE);
  saveRDS(std[order(rownames(std)),], file=paste(path, 'r/study.rds', sep='/')); 
  
  ## full analysis info by ID
  a<-cbind(ana[rownames(an), ], an, PubMed=sub('pmid', '', as.vector(ana[rownames(an), 'Study_ID'])));
  cnm<-c( # fields to be included
    Name = 'Name',
    'Study ID' = 'Study_ID',
    URL = 'URL',
    'PubMed' = 'PubMed',
    'Number of variants' = 'Num_Variants',
    'Minimum Phred score' = 'Min_Phred',
    'Maximum Phred score' = 'Max_Phred',
    Trait = 'DISEASE_TRAIT',
    Journal = 'JOURNAL',
    Date = 'DATE',
    Author = 'FIRST_AUTHOR',
    Platform = 'PLATFORM_SNPS_PASSING_QC',
    'Initial sample' = 'INITIAL_SAMPLE_SIZE',
    'Replication sample' = 'REPLICATION_SAMPLE_SIZE'
  )
  cnm<-cnm[cnm %in% colnames(a)];
  a<-a[, cnm];
  for (i in 1:ncol(a)) a[[i]]<-as.vector(a[[i]]);
  a[is.na(a)]<-'';
  colnames(a)<-names(cnm);
  # create list
  id2ana<-lapply(rownames(a), function(id) {
    flds<-a[id, ];
    sid<-as.vector(ana[id, 'Study_ID']);
    as.list(c(ID = id, Source='GWAScatalog', flds, 'Study Name' = as.vector(std[sid, 'Name']), 'Study full description' = as.vector(std[sid, 'Description'])));
  });
  names(id2ana)<-rownames(a);
  id2ana<-id2ana[order(names(id2ana))];
  saveRDS(id2ana, file=paste(path, 'r/analysis_by_id.rds', sep='/'))
  
  ## full study info by ID
  id2std<-lapply(rownames(std), function(id) {
    list(
      ID = id,
      Source = 'GWAScatalog',
      Name = std[id, 'Name'],
      URL = std[id, 'URL'],
      PubMed = sub('pmid', '', id),
      'Number of analyses' = std[id, 'Num_Analyses'],
      Analyses = rownames(ana)[ana$Study_ID==id & !is.na(ana$Study_ID)],
      Journal = pubmed[sub('^pmid', '', id), 'Journal'],
      Date = pubmed[sub('^pmid', '', id), 'Date'],
      'Full description' = std[id, 'Description']
    )
  });
  names(id2std)<-rownames(std);
  saveRDS(id2std[order(names(id2std))], file=paste(path, 'r/study_by_id.rds', sep='/'));
  
  ana2pm<-lapply(id2ana, function(x) x$PubMed);
  std2pm<-lapply(id2std, function(x) x$PubMed);
  saveRDS(ana2pm, file=paste(path, 'r/analysis2pubmed.rds', sep='/'));
  saveRDS(std2pm, file=paste(path, 'r/study2pubmed.rds', sep='/'));
  
  # Mapping to genes
  gw1 <- gw0[, c("PUBMEDID", "DISEASE_TRAIT", "UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID", "SNP_GENE_IDS", 
             "UPSTREAM_GENE_DISTANCE", "DOWNSTREAM_GENE_DISTANCE", "SNP_ID_CURRENT")];
  gw1 <- gw1[gw1$UPSTREAM_GENE_DISTANCE<=10^6 | is.na(gw1$UPSTREAM_GENE_DISTANCE), ];
  gw1 <- gw1[gw1$DOWNSTREAM_GENE_DISTANCE<=10^6 | is.na(gw1$DOWNSTREAM_GENE_DISTANCE), ];
  dis <- gsub(' ', '_', as.vector(gw1$DISEASE_TRAIT)); 
  dis <- gsub('^_', '', dis);
  dis <- gsub('_$', '', dis); 
  gns <- gw1[, c("UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID", "SNP_GENE_IDS")];
  gns <- paste(gns[, 1], gns[, 2], gns[, 3], sep=', ');
  gns <- strsplit(gns, ', ');
  gns <- lapply(gns, function(g) g[g!='']); 
  mp1 <- split(gns, gw1$PUBMEDID);
  mp1 <- lapply(mp1, function(x) unique(unlist(x))); 
  mp1 <- lapply(mp1, function(x) x[x!='' & !is.na(x)]); 
  mp1 <- mp1[sapply(mp1, length)>0];
  mp2 <- split(gns, dis); 
  mp2 <- lapply(mp2, function(x) unique(unlist(x))); 
  mp2 <- lapply(mp2, function(x) x[x!='' & !is.na(x)]); 
  mp2 <- mp2[sapply(mp2, length)>0];  
  mp3 <- split(gns, paste(gw1$PUBMEDID, dis, sep='_'));
  mp3 <- lapply(mp3, function(x) unique(unlist(x))); 
  mp3 <- lapply(mp3, function(x) x[x!='' & !is.na(x)]); 
  mp3 <- mp1[sapply(mp3, length)>0];  
  saveRDS(mp1, file=paste(path, 'r/pubmed2gene.rds', sep='/'));
  saveRDS(mp2, file=paste(path, 'r/disease2gene.rds', sep='/'));
  saveRDS(mp3, file=paste(path, 'r/pubmed_disease2gene.rds', sep='/'));
  
  # Mapping to SNPs
  snp <- as.vector(gw0$SNP_ID_CURRENT);
  snp <- paste('rs', snp, sep=''); 
  snp[grep('x', snp)] <- '';
  dis <- gsub(' ', '_', as.vector(gw0$DISEASE_TRAIT)); 
  dis <- gsub('^_', '', dis);
  dis <- gsub('_$', '', dis); 
  mp1 <- split(snp, gw0$PUBMEDID);
  mp1 <- lapply(mp1, function(x) unique(unlist(x))); 
  mp1 <- lapply(mp1, function(x) x[x!='' & !is.na(x)]); 
  mp1 <- mp1[sapply(mp1, length)>0];
  mp2 <- split(snp, dis); 
  mp2 <- lapply(mp2, function(x) unique(unlist(x))); 
  mp2 <- lapply(mp2, function(x) x[x!='' & !is.na(x)]); 
  mp2 <- mp2[sapply(mp2, length)>0];  
  mp3 <- split(snp, paste(gw0$PUBMEDID, dis, sep='_'));
  mp3 <- lapply(mp3, function(x) unique(unlist(x))); 
  mp3 <- lapply(mp3, function(x) x[x!='' & !is.na(x)]); 
  mp3 <- mp1[sapply(mp3, length)>0];  
  saveRDS(mp1, file=paste(path, 'r/pubmed2snp.rds', sep='/'));
  saveRDS(mp2, file=paste(path, 'r/disease2snp.rds', sep='/'));
  saveRDS(mp3, file=paste(path, 'r/pubmed_disease2snp.rds', sep='/'));
  
  list(
    snp = rownames(tbl),
    analysis = rownames(ana),
    study = rownames(std),
    pubmed = rownames(pubmed)
  )
}