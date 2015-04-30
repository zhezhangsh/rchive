# Download and process GWAS results from GWAS catalog
ParseGwasCatalog<-function(url, path=paste(RCHIVE_HOME, 'data/gwas/public/gwascatalog', sep='/'), update.all.pubmed=FALSE) {
  # url    Source file URL
  # path   Path to output files
  # update.all.pubmed   Re-download all PubMed entries if TRUE
  
  library(RCurl);
  library(NCBI2R);
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  # download the table 
  ln<-strsplit(getURL(url), '\n')[[1]];
  gw<-do.call('rbind', strsplit(ln, '\t'));
  colnames(gw)<-gw[1,];
  gw<-gw[-1, ];
  gw<-data.frame(gw, stringsAsFactors=FALSE);
  
  # this will download an outdated version of GWAScatalog, 
  # but its column names were used across the rest of this function
  # so it's downloaded just to provide the column names
  gw0<-GetPublishedGWAS(); 
  if (ncol(gw) == ncol(gw0)) colnames(gw)<-colnames(gw0) else warning("Unmatched column names between versions!")
  
  # save a copy of full table
  fn.src<-sapply(strsplit(url, '://'), function(x) x[2]);
  fn.src<-gsub('/', '_', fn.src);
  saveRDS(gw, file=paste(fn.src, '.rds', sep=''));
  
  ##########################################################################################################
  # start working on the table
  
  # fields of numeric values
  numeric.fields<-c('Chr_id', 'Chr_pos', 'Merged', 'Intergenic', 'RiskAlleleFrequency', 'pValue', 'Pvalue_mlog', 'ORorbeta'); # fields shold be numeric
  num<-gw[, numeric.fields];
  num<-apply(num, 2, as.numeric);
  gw[, numeric.fields]<-num;
  
  # fields must have valid value
  required.fields<-c('Snp_id_current', 'Chr_id', 'Chr_pos', 'Pvalue_mlog');
  req<-gw[, required.fields];
  n<-apply(req, 2, function(x) is.na(x) | x=='');
  cnv<-gw[gw$CNV=='Y', , drop=FALSE]; # CNV variants
  inc<-gw[rowSums(n)>0, , drop=FALSE];
  gw<-gw[gw$CNV!='Y' & rowSums(n)==0, , drop=FALSE]; 
  saveRDS(cnv, file=paste(path, 'r', 'gwascatalog_just_cnv.rds', sep='/'));
  saveRDS(inc, file=paste(path, 'r', 'gwascatalog_incomplete_info.rds', sep='/'));
  saveRDS(gw, file=paste(path, 'r', 'gwascatalog_complete_info.rds', sep='/'));
  
  #####################
  gw<-data.frame(id=paste('rs', gw[['Snp_id_current']], sep=''), phred=as.integer(round(10*abs(gw[['Pvalue_mlog']]))), gw, stringsAsFactors=FALSE);
  #####################
  
  # rename chromosomes
  chr<-as.character(gw[['Chr_id']]);
  chr[chr=='23']<-'X';
  chr[chr=='24']<-'Y';
  chr[chr=='25']<-'MT';
  gw[['Chr_id']]<-chr;
  
  # the combination of SNP ID, PubMed ID, Trait, and p value text should be unique
  dup.fields<-c('ReportedGenes', 'Mapped_gene', 'Region', 'Snp_gene_ids'); # the fields that might cause duplicated entries
  gw<-gw[!duplicated(gw[, !(colnames(gw) %in% dup.fields)]), ];
  uniq<-paste(gw[['id']], gw[['PUBMEDID']], gw[['DiseaseTrait']], gw[['pValuetext']], sep='@_@');
  dup<-uniq[duplicated(uniq)];# duplicated entries
  dups<-lapply(dup, function(dup) gw[uniq==dup, ]);
  dups<-sapply(dups, function(d) rownames(d[order(d[, 'phred'])])[nrow(d)]); # pick one with the highest phred score
  gw<-rbind(gw[!(uniq %in% dup), ], gw[dups,]); # duplicated enries removed
  gw<-gw[!is.na(gw$phred), ]; # Entry will be removed if reported p value is 0
  
  ###############################################################################################################
  ###############################################################################################################
  cat('Re-formatting data ...\n');
  
  # re-format p value text
  txt<-gw[['pValuetext']];
  txt[is.na(txt)]<-'';
  txt<-gsub(' ', '', txt);
  txt<-gsub("[()]", '', txt);
  #anno[,4]<-txt;
  
  # unique analysis and study
  ana<-paste(gw[['PUBMEDID']], gw[['DiseaseTrait']], gw[['pValuetext']], sep='@_@'); # analyses
  sty<-split(ana, gw[['PUBMEDID']]); # studies
  sty<-lapply(sty, unique);
  sty<-lapply(sty, sort);
  
  # analysis ID
  ana.nm<-unlist(sty, use.names=FALSE); # analysis names
  ana.id<-paste('pmid', rep(names(sty), sapply(sty, length)), '_', unlist(sapply(sty, function(x) 1:length(x))), sep=''); # analysis ID
  names(ana.id)<-ana.nm;
  id<-as.vector(ana.id[ana]); # analysis ID of each row
  
  # SNP positions
  pos<-gw[, c('id', 'Chr_id', 'Chr_pos')];
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
  #an<-gw[!duplicated(id), c(4, 7, 6, 8, 35, 9, 10, 11, 12)]; # select fields
  an<-gw[!duplicated(id), ]; # one analysis, one row
  rownames(an)<-id[!duplicated(id)];
  sty.id<-paste('pmid', an$PUBMEDID, sep=''); # Non-unique study ID
  st.id<-unique(sty.id); # unique study ID
  
  # Get article info from pubmed
  if (!update.all.pubmed & file.exists(paste(path, 'r/pubmed_downloaded.rds', sep='/'))) { # existing info
    print('Loading existing pubmed information ...');
    pm<-readRDS(paste(path, 'r/pubmed_downloaded.rds', sep='/'));
    id<-sub('^pmid', '', st.id);
    new.pub<-id[!(id %in% names(pm))];
    if (length(new.pub) > 0) { # there are new pubmed articles
      new.pm<-lapply(new.pub, GetPubMed);
      names(new.pm)<-new.pub;
      saveRDS(new.pm, file=paste(path, 'r/pubmed_new_article.rds', sep='/'));
      pm<-c(pm, new.pm);
      saveRDS(pm,file=paste(path, 'r/pubmed_downloaded.rds', sep='/'));
    }
  } else {
    library(NCBI2R)
    pm<-lapply(sub('^pmid', '', st.id), GetPubMed);
    names(pm)<-sub('^pmid', '', st.id);
    saveRDS(pm,file=paste(path, 'r/pubmed_downloaded.rds', sep='/'));
  }    
  
  # PubMed table
  pubmed<-sapply(pm, function(pm) pm[1, c('TI', 'JT', 'DA', 'link', 'AB')]);
  pubmed<-data.frame(t(pubmed), stringsAsFactors=FALSE, row.names=names(pm));
  colnames(pubmed)<-c('Title', 'Journal', 'Date', 'URL', 'Abstract');
  for (i in 1:ncol(pubmed)) pubmed[[i]]<-as.vector(unlist(pubmed[[i]]));
  saveRDS(pubmed, file=paste(path, 'r/pubmed.rds', sep='/'));
  
  # PubMed list named by PMID
  pm.by.id<-lapply(pm, function(pm) as.list(pm[1,]));
  names(pm.by.id)<-names(pm);
  saveRDS(pm.by.id, file=paste(path, 'r/pubmed_by_id.rds', sep='/'));
  
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
  names(pm)<-paste('pmid', names(pm), sep='');
  nm<-as.vector(sapply(pm, function(pm) pm$TI));
  url<-as.vector(sapply(pm, function(pm) pm$link));
  desc<-as.vector(sapply(pm, function(pm) pm$AB));
  n.ana<-as.vector(as.vector(table(ana$Study)[names(pm)]));
  n.ana[is.na(n.ana)]<-0;
  std<-data.frame(Source=rep('GWAScatalog', length(pm)), Name=nm, Num_Analyses=n.ana, URL=url, Description=desc, row.names=names(pm));
  std<-std[order(rownames(std)),];
  names(pm)<-sub('^pmid', '', names(pm));
  saveRDS(std, file=paste(path, 'r/study.rds', sep='/'))
  
  ## full analysis info by ID
  a<-cbind(ana[rownames(an), ], an);
  cnm<-c( # fields to be included
    Name = 'Name',
    'Study ID' = 'Study_ID',
    URL = 'URL',
    'Number of variants' = 'Num_Variants',
    'Minimum Phred score' = 'Min_Phred',
    'Maximum Phred score' = 'Max_Phred',
    Trait = 'DiseaseTrait',
    Journal = 'Journal',
    Date = 'Date',
    Author = 'FirstAuthor',
    Platform = 'Platform',
    'Initial sample' = 'InitialSampleDescription',
    'Replication sample' = 'ReplicationSampleDescription'
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
  sid<-as.vector(ana[rownames(an), ]$Study_ID);
  jnl<-sapply(split(as.vector(an$Journal), sid), function(j) paste(unique(j), collapse=';'));
  dt<-sapply(split(as.vector(an$Date), sid), function(j) paste(unique(j), collapse=';'));
  for (i in 1:ncol(std)) std[[i]]<-as.vector(std[[i]]);
  id2std<-lapply(rownames(std), function(id) {
    list(
      ID = id,
      Source = 'GWAScatalog',
      Name = std[id, 'Name'],
      URL = std[id, 'URL'],
      'Number of analyses' = std[id, 'Num_Analyses'],
      Journal = as.vector(jnl[id]),
      Date = as.vector(dt[id]),        
      'Full description' = std[id, 'Description']
    )
  });
  names(id2std)<-rownames(std);
  id2std<-id2std[order(names(id2std))];
  saveRDS(id2std, file=paste(path, 'r/study_by_id.rds', sep='/'));
  
  list(
    snp = rownames(tbl),
    analysis = rownames(ana),
    study = rownames(std),
    pubmed = rownames(pubmed)
  )
}