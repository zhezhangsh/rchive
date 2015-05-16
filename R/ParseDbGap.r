# Download and process GWAS results from dbGap
DownloadDbGap<-function(url="ftp://ftp.ncbi.nlm.nih.gov/dbgap//Analysis_Table_of_Contents.txt", 
                     path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'), 
                     update.all.analysis=FALSE,
                     redundant.studies=c('216', '88', '342')) {
  # url                   File with analysis metadata and download location
  # path                  Path to output files
  # update.all.analysis   Re-download all PubMed entries if TRUE, takes very long time to do so
  # redundant.studies     Exclude redundant studies of the same analyses belong to other studys
  
  library(RCurl);
  library(NCBI2R);
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if (!file.exists(paste(path, 'r/full_table', sep='/'))) dir.create(paste(path, 'r/full_table', sep='/'), recursive=TRUE);
  if (!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
    
  # Original metadata with each ftp file of analysis; parsed into a table
  fn.meta<-paste(path, 'r', 'original_metadata.rds', sep='/');
  if (file.exists(fn.meta)) meta<-readRDS(fn.meta) else {
    emp<-character(0);
    meta<-data.frame(name=emp, study=emp, genome=emp, dbsnp=emp, description=emp, method=emp, file=emp, stringsAsFactors=FALSE);
    #saveRDS(meta, file=fn.meta);
  }
  # Original metadata with each ftp file of analysis; as a list
  fn.meta.lst<-paste(path, 'r', 'original_metadata_list.rds', sep='/');
  if (file.exists(fn.meta.lst)) meta.lst<-readRDS(fn.meta.lst) else {
    meta.lst<-list();
    #saveRDS(meta.lst, file=fn.meta.lst);
  }
  
  # Download and parse metadata table
  ln<-strsplit(getURL(url), '\n')[[1]];
  gw<-do.call('rbind', strsplit(ln, '\t'));
  colnames(gw)<-gw[1,];
  gw<-gw[-1, , drop=FALSE];
  gw<-data.frame(gw, stringsAsFactors=FALSE);
  gw<-gw[!(gw[,1] %in% redundant.studies), ];
  
  ######################################################################################
  # Download results of individual analyses
  # file information
  fn.ftp<-gw$analysis.result;
  fn.gz<-sapply(strsplit(fn.ftp, '/'), function(x) x[length(x)]);
  id.study<-sapply(strsplit(fn.gz, '\\.'), function(x) x[1]);
  id.analysis<-sapply(strsplit(fn.gz, '\\.'), function(x) x[2]);
  fn.gz<-paste(path, 'src', fn.gz, sep='/');
  fn.r<-paste(path, 'r/full_table', paste(id.analysis, '.rds', sep=''), sep='/');
  file.info<-cbind(analysis=id.analysis, study=id.study, ftp=fn.ftp, gz=fn.gz, r=fn.r);
  
  ####################################################################################
  # load in table and save to local
  loaded<-sub('.rds', '', sapply(strsplit(fn.r, '/'), function(x) x[length(x)]))
  no.meta<-!(loaded %in% rownames(meta)) | !(loaded %in% names(meta.lst)); # missing metadata info.
  if (update.all.analysis) no.tbl<-file.info else no.tbl<-file.info[!file.exists(file.info[, 'r']) | no.meta, , drop=FALSE];
  if (nrow(no.tbl) > 0) {
    hd<-lapply(1:nrow(no.tbl), function(i) {
      x<-no.tbl[i,];
      cat("loading analysis", x[1], '\n');
      if (!file.exists(x[4])) download.file(x[3], x[4]); 
      ln<-scan(x[4], what='', flush=TRUE, sep='\n');
      
      # very tediously parse file header and add header info to metadata table and list
      hd<-ln[substr(ln, 1, 2)=='# '];
      hd<-sub('^# ', '', hd);
      hd<-strsplit(hd, ':\t');
      hd<-lapply(hd, function(hd) sub('^ ', '', hd));
      hd<-lapply(hd, function(hd) sub(' $', '', hd));      
      hd<-do.call('rbind', hd);
      hd.nm<-tolower(hd[,1]);
      hd<-hd[,2];
      names(hd)<-hd.nm;
      
      # save full result table
      cnm<-strsplit(ln[length(hd)+1], '\t')[[1]];
      tbl<-strsplit(ln[(length(hd)+2):length(ln)], '\t');
      tbl<-t(sapply(tbl, function(x) x[1:length(cnm)]));
      tbl<-data.frame(tbl, stringsAsFactors=FALSE); # all columns are character
      colnames(tbl)<-cnm;
      saveRDS(tbl, x[5]);
      
      hd;
    })
    
    meta.lst[no.tbl[,1]]<-hd;
    new.rows<-sapply(1:length(hd), function(i) {
      hd<-hd[[i]];
      as.vector(c(hd['name'], no.tbl[i, 2], hd[c('human genome build', 'dbsnp build', 'description', 'method')], no.tbl[i, 3]));
    });
    new.rows<-t(new.rows);
    rownames(new.rows)<-no.tbl[,1];
    new.rows<-new.rows[!(rownames(new.rows) %in% rownames(meta)), , drop=FALSE];
    if (nrow(new.rows) > 0) {
      df<-data.frame(new.rows, stringsAsFactors=FALSE);
      colnames(df)<-colnames(meta);
      meta<-rbind(meta, df);
    }
  } 
  ####################################################################################
  
  saveRDS(meta, file=fn.meta);
  saveRDS(meta.lst, file=fn.meta.lst);
  
  meta;
}

# Retrieve a specified test statistics: p value, effect size, or allele frequency from downloaded dbGaP results files
RetrieveDbGapStat<-function(id.ana, id.std, stat.name=c('p value', 'effect size', 'allele frequency'), 
                            path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'),
                            own.table.min=20) {
  # id.ana, id.std        Matched analysis and study IDs, preferentially including full list of dbGaP analysis, both previously loaded and unloaded analysis
  # stat.name             The type of test statistics to summarize. Integer (1, 2, 3) or name of the statistics
  # path                  Path to input/output files  
  # own.table.min         Minimum number of analyses to make a study its own table
  
  library(rchive);
  
  own.table.min<-max(3, own.table.min); # At lease 3 analyses to make their own table
  
  # specify the type and name of test statistics to summarize
  stat.types<-c('p value' = 'phred', 'effect size' = 'es', 'allele frequency' = 'maf');
  stat.type<-stat.types[stat.name[1]];
  if (is.na(stat.type)) stat.type<-c('p value'=1);
  
  # Column names that fit to the type of test statistics
  cnm<-list(
    'phred' = c('p-value', 'p_value', 'p value', 'pvalue', 'p'),
    'es' = c('effect size', 'effect_size', 'effect-size', 'or', 'odds ratio', 'odds-ratio', 'odds_ratio', 'beta', 'gee', 'gee beta', 'gee-beta', 'gee_beta'),
    'maf' = c('maf', 'af', 'minor allele frequency', 'minor_allele_frequency', 'minor-allele-frequency', 'allele frequency', 'allele_frequency', 'allele-frequency')
  )[stat.type];
  
  # Group analysis by study IDs
  if (length(id.ana) != length(id.std)) stop("Number of analyses and number of studies don't match!");
  names(id.ana)<-id.std;
  id.ana<-id.ana[!is.na(id.ana) & !is.na(id.std) & id.ana!='' & id.std!=''];
  std2ana<-split(as.vector(id.ana), id.std);
  std2ana<-lapply(std2ana, function(ids) {
    fn<-paste(path, '/r/full_table/', ids, '.rds', sep='');
    names(fn)<-ids;
    fn[file.exists(fn)];
  });
  std2ana<-std2ana[sapply(std2ana, length)>0];
  
  if (length(std2ana)>0) {
    # the table save the analyses of studies don't have their own table
    fn.tbl0<-paste(path, '/r/', names(cnm)[1], '_', 'phs000000.rds', sep='');
    if (file.exists(fn.tbl0)) tbl0<-readRDS(fn.tbl0) else tbl0<-matrix(nr=0, nc=0);
    ana0<-colnames(tbl0);
    
    fn.std<-paste(path, '/r/', names(cnm)[1], '_', names(std2ana), '.rds', sep=''); # file name of the tables storing stats of a study
    names(fn.std)<-names(std2ana);
    existing<-file.exists(fn.std); # whether the table corresponding to the study already exists
    std2ana.existing<-std2ana[existing]; # studies already have a table 
    std2ana.new<-std2ana[!existing]; # studies have no existing table 

    # add new analysis to existing tables
    if (length(std2ana.existing) > 0) {
      sapply(names(std2ana.existing), function(nm) {
        cat('Updating study', nm, '\n');
        tbl<-readRDS(fn.std[nm]);
        ana<-std2ana.existing[[nm]];
        ana<-ana[!(names(ana) %in% colnames(tbl))];
        if (length(ana) > 0) {
          tbl1<-LoadColumn(ana, cnm[[1]], row=c('snp id', 'snp_id', 'snp-id'), collapse='min', is.numeric=TRUE);
          tbl<-MergeTable(list(tbl, tbl1));
          saveRDS(tbl[!is.na(rowMeans(tbl, na.rm=TRUE)), !is.na(colMeans(tbl, na.rm=TRUE)), drop=FALSE], file=fn.std[nm]);
        }
      });
    }
    
    # new analysis to new tables
    if (length(std2ana.new) > 0) {
      tbl<-lapply(names(std2ana.new), function(nm) {
        cat('Updating study', nm, '\n');
        ana<-std2ana.new[[nm]];
        tbl.loaded<-tbl0[, colnames(tbl0) %in% names(ana), drop=FALSE];
        ana.unloaded<-ana[!(names(ana) %in% colnames(tbl0))];
        if (length(ana.unloaded) > 0) {
          tbl<-LoadColumn(std2ana.new[[nm]], cnm[[1]], row=c('snp id', 'snp_id', 'snp-id'), collapse='min', is.numeric=TRUE);
          tbl<-round(-10*log10(tbl));
          tbl<-MergeTable(list(tbl, tbl.loaded));
        } else tbl<-tbl.loaded;
        tbl;
      });
      
      # Save to new file or merge to table for orphan analyses
      for (i in 1:length(tbl)) {
        if (ncol(tbl[[i]]) >= own.table.min) {
          t<-tbl[[i]]; 
          t<-t[!is.na(rowMeans(t, na.rm=TRUE)), !is.na(colMeans(t, na.rm=TRUE)), drop=FALSE];
          saveRDS(t, fn.std[names(std2ana.new)[i]]);
        } else tbl0<-MergeTable(list(tbl0, tbl[[i]]));
      }
      
      # Remove analyses having their own table from orphan table
      not.orphan<-unlist(lapply(tbl, function(t) if (ncol(t) > own.table.min) colnames(t) else c()), use.names=FALSE);
      not.orphan<-unique(c(not.orphan, unlist(lapply(std2ana.existing, names), use.names=FALSE)));
      tbl0<-tbl0[, !(colnames(tbl0) %in% not.orphan), drop=FALSE];
      
      if (!setequal(ana0, colnames(tbl0))) saveRDS(tbl0[!is.na(rowMeans(tbl0, na.rm=TRUE)), !is.na(colMeans(tbl0, na.rm=TRUE)), drop=FALSE], file=fn.tbl0);
    }
  }
}

################################################################################
# get HTML page of the study from dbGAP, and retrieve fields
GetDbGapStudy<-function(id) {
  #id     dbGaP study ID
  
  # Get the page of the first version of the study
  url<-paste('http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=', id, '.v1.p1', sep='');
  download.file(url, id);
  html<-scan(id, sep='\n', what='', flush=TRUE);
  if (length(html[html=="The requested study does not exist in the dbGaP system."])>0) { # original version not available
    url<-paste('http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=', id, '.v1.p2', sep='');
    download.file(url, id);
    html<-scan(id, sep='\n', what='', flush=TRUE);
  }
  
  
  # check if there is newer version
  l<-html[grep('^<a', html)];
  l<-l[grep(paste('?study_id=', id, '.v', sep=''), l)];
  
  # replace page with newer version if there is one
  if (length(l)>0) {
    l<-sort(l)[length(l)];
    l<-strsplit(l, '[\"=]')[[1]];
    v<-l[grep(paste('^', id, '.v', sep=''), l)];
    url<-paste('http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=', v, sep='');
    download.file(url, id);
    html<-scan(id, sep='\n', what='', flush=TRUE);
  }
  
  ############### Get Pubmed ID associated to this study
  #ln<-html[grep('>Pubmed<', html)];
  #ln<-ln[grep('^<a href=', ln)][1];
  #if (length(ln) == 0) pmid<-NA else {
  #    url<-strsplit(ln, '"')[[1]][[2]];
  
  #}
  
  diseases<-html[grep('DB=mesh', html)+1];
  diseases<-paste(unique(diseases), collapse='; ');
  
  # Get title
  html<-html[(grep("<span id=\"study-name\" name=\"study-name\">", html)+1):length(html)];
  title<-html[1];
  
  # remove HTML tags etc.
  desc<-gsub('\t', ' ', html);
  desc<-gsub("<.*?>", "", desc); # remove html tags
  desc<-gsub("^ *|(?<= ) | *$", "", desc, perl = TRUE); # remove extra spaces.
  desc<-gsub('^[ \t|]', '', desc); 
  desc<-gsub('[ \t|]$', '', desc);
  desc<-desc[desc!='']; 
  desc<-desc[!(desc %in% c("Important Links and Information", "Request access via Authorized Access", "Instructions for requestors", "Data Use Certification (DUC) Agreement", "Acknowledgement Statements", "Participant Protection Policy FAQ"))];
  ty<-desc[(grep('Study Type:', desc)[1]+1)];
  hs<-desc[(grep('Study History', desc)[1]+1)];
  pi<-desc[(grep('Principal Investigator', desc)[1]+1)];
  dsc<-desc[(grep('Study Description', desc)[1]+1)];
  
  file.remove(id);
  out<-list(title=title, url=url, disease=diseases, type=ty, history=hs, pi=pi, description=dsc);
  as.list(sapply(out, function(out) if (length(out)==0) '' else out));
}
################################################################################


# Summarize metadata infomration of dbGaP studies
SummarizeDbGap<-function(meta, path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'), 
                         update.all.study=FALSE, update.all.pubmed=FALSE) {
  # meta              Analysis metadata table prepared by the DownloadDbGap function
  # path              Path to output folder
  # update.all.study  If TRUE, re-download information of all studies from source; Just download the new ones if FASLE
  # update.all.pubmed If TRUE, update all pubmed-related information; otherwise, using previously saved information when it's available
    
  std2ana<-split(rownames(meta), meta$study);
    
  ###############################################
  # Basic phred statistics of analysis table 
  fn.stat<-paste(path, 'r/analysis_stat.rds', sep='/'); # save analysis basic stats
  if (file.exists(fn.stat)) {
    stat<-readRDS(fn.stat);
    if (setequal(rownames(stat), unlist(std2ana, use.names=FALSE))) do.stat<-FALSE else do.stat<-TRUE;
  } else do.stat<-TRUE;
  
  if (do.stat) {
    f<-dir(paste(path, 'r', sep='/'));
    f<-f[grep('^phred', f)];
    f<-f[grep('.rds$', f)];
    f<-paste(path, 'r', f, sep='/');
    stat<-lapply(f, function(f) {
      print(f);
      mtrx<-readRDS(f);
      apply(mtrx, 2, function(p) {
        c(length(p[!is.na(p)]), min(p, na.rm=TRUE), max(p, na.rm=TRUE))
      })
    });
    stat<-t(do.call('cbind', stat));
    colnames(stat)<-c('N', 'Min', 'Max');
    saveRDS(stat, file=fn.stat);
  }
  
  std2ana<-lapply(std2ana, function(x) sort(intersect(x, rownames(stat))));
  std2ana<-std2ana[sapply(std2ana, length)>0];
  if (length(std2ana) == 0) stop("No valid study");

  # Full study annotation table
  fn.dbgap<-paste(path, 'r/study_downloaded.rds', sep='/'); # save annotation of studies retrievd from dbGap
  st.id<-names(std2ana);
  if (file.exists(fn.dbgap) & !update.all.study) {
    st<-readRDS(fn.dbgap);
    st<-st[names(st) %in% names(std2ana)];
    new.id<-setdiff(st.id, names(st));
    if (length(new.id) > 0) {
      st0<-lapply(new.id, GetDbGapStudy);
      names(st0)<-new.id;
      st<-c(st, st0);
      st<-st[order(names(st))];
      saveRDS(st0, file=sub('downloaded.rds', 'new.rds',fn.dbgap));
      saveRDS(st, file=fn.dbgap);
    }
  } else {
    st<-lapply(st.id, function(x) {print(x); GetDbGapStudy(x); }); # get study info from dbGAP web page
    names(st)<-st.id;
    st<-st[order(names(st))];
    saveRDS(st, file=fn.dbgap);
  }
  
  ###############################################
  # study anotation table
  nm<-as.vector(sapply(st, function(st) st$title));
  url<-as.vector(sapply(st, function(st) st$url));
  desc<-as.vector(sapply(st, function(st) st$description));
  n.ana<-sapply(std2ana, length);
  n.ana[is.na(n.ana)]<-0;
  std<-data.frame(Source=rep('dbGaP', length(st)), Name=nm, Num_Analyses=n.ana, URL=url, Description=desc, row.names=names(st));
  std<-std[order(rownames(std)),];
  saveRDS(std, file=paste(path, 'r/study.rds', sep='/'))
  
  # Add summary statistics to analysis table
  ids<-intersect(rownames(meta), rownames(stat));
  an<-meta[ids, , drop=FALSE];
  url0<-as.vector(std$URL); # build the analysis URL
  names(url0)<-rownames(std);
  url<-paste(url0[as.vector(an$study)], '&pha=', as.numeric(sub('pha', '', rownames(an))), sep='');
  url<-sub('study.cgi', 'analysis.cgi', url);
  ana<-data.frame(Source=rep('dbGaP', nrow(an)), Name=as.vector(an$name), Study_ID=as.vector(an$study), URL=url,
                  Num_Variants=stat[, 'N'][ids], Min_Phred=stat[, 'Min'][ids], Max_Phred=stat[, 'Max'][ids], 
                  Description=an$description, row.names=rownames(an));
  ana<-ana[order(rownames(ana)), ];
  saveRDS(ana, file=paste(path, 'r/analysis.rds', sep='/'))
  
  ###############################################
  ## full study info by ID
  sid<-as.vector(ana[rownames(an), ]$Study_ID);
  for (i in 1:ncol(std)) std[[i]]<-as.vector(std[[i]]);
  id2std<-lapply(rownames(std), function(id) {
    list(
      ID = id,
      Source = 'dbGaP',
      Name = std[id, 'Name'],
      URL = std[id, 'URL'],
      'Number of analyses' = std[id, 'Num_Analyses'],
      'Principal investigator' = st[[id]]$pi,
      'Diseases' = st[[id]]$disease,
      'Study type' = st[[id]]$type,
      'Study history' = st[[id]]$history,
      'Full description' = std[id, 'Description']
    )
  });
  names(id2std)<-rownames(std);
  id2std<-id2std[order(names(id2std))];
  saveRDS(id2std, file=paste(path, 'r/study_by_id.rds', sep='/'))
  
  ###############################################
  ## full analysis info by ID
  a<-cbind(ana[rownames(an), ], an);
  a$file<-paste('ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies', a$file, sep='/');
  
  cnm<-c( # fields to be included
    Name = 'Name',
    'Study ID' = 'Study_ID',
    URL = 'URL',
    'Number of variants' = 'Num_Variants',
    'Minimum Phred score' = 'Min_Phred',
    'Maximum Phred score' = 'Max_Phred',
    'Genome version' = 'genome',
    'dbSNP version' = 'dbsnp',
    'Source file' = 'file',
    'Method' = 'method',
    'Description' = 'description'
  );
  cnm<-cnm[cnm %in% colnames(a)];
  a<-a[, cnm];
  for (i in 1:ncol(a)) a[[i]]<-as.vector(a[[i]]);
  a[is.na(a)]<-'';
  colnames(a)<-names(cnm);
  # create list
  id2ana<-lapply(rownames(a), function(id) {
    flds<-as.list(a[id, ]);
    sid<-as.vector(ana[id, 'Study_ID']);
    c(ID = id, Source='dbGap', flds);
  });
  names(id2ana)<-rownames(a);
  id2ana<-id2ana[order(names(id2ana))];
  saveRDS(id2ana, file=paste(path, 'r/analysis_by_id.rds', sep='/'));
  
  ###################################################################################################################
  ###################################################################################################################
  # Now, associate analysis and study to Pubmed articles
  std<-readRDS(file=paste(path, 'r/study.rds', sep='/'));
  ana<-readRDS(file=paste(path, 'r/analysis.rds', sep='/'));
  id2std<-readRDS(file=paste(path, 'r/study_by_id.rds', sep='/'));
  id2ana<-readRDS(file=paste(path, 'r/analysis_by_id.rds', sep='/'));
  
  library(RCurl);
  
  ##################################################################
  # get one page of Pubmed IDs based on given URL
  getPMID<-function(url) {
    url<-url[1];
    print(url);
    if (url=='' | is.na(url)) '' else {
      html<-strsplit(getURL(url), '\n')[[1]];
      ln<-html[grep('^<a href=', html)];
      ln<-ln[grep('db=pubmed', ln)];
      ln<-strsplit(ln, '[=\\[]');
      id<-sapply(ln, function(ln) ln[length(ln)-1]);
      id;
    }
  }
  ##################################################################
  
  
  ############
  # Analysis #
  ############
  
  ###############################################
  url0<-as.vector(ana$URL); # URL to each analysis
  names(url0)<-rownames(ana);
  
  # Get URL to CGI script that reports pubmed linked to each analysis
  if (!update.all.pubmed & file.exists(paste(path, 'r/url_analysis2pubmed.rds', sep='/'))) {
    url1<-readRDS(file=paste(path, 'r/url_analysis2pubmed.rds', sep='/'));
    url0<-url0[!(names(url0) %in% names(url1))];
  } else {
    url1<-c();
  }
  if (length(url0) > 0) {
    url.new<-sapply(url0, function(url) {
      html<-strsplit(getURL(url), '\n')[[1]];
      ln<-html[grep('initializeAnalysisReferences', html)][1];
      ln<-strsplit(ln, '[\";]')[[1]];
      ln<-ln[grep('^initializeAnalysisReferences', ln)][1]
      ln<-strsplit(ln, '[(\',)]')[[1]]; 
      print(ln[3]);   
      if (length(ln) !=7) '' else {
        paste("http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetAnalysisReference.cgi?pha=", ln[3], '&version=', ln[5], '&page_number=1', sep='');
      }
    });
    names(url.new)<-names(url0);
    url1<-c(url1, url.new);
  }
  url1<-url1[rownames(ana)];
  saveRDS(url1, file=paste(path, 'r/url_analysis2pubmed.rds', sep='/'));
  
  ###############################################
  # Get Pubmed mapped to the analysis from a cgi script of dbGap
  url1<-readRDS(file=paste(path, 'r/url_analysis2pubmed.rds', sep='/'));
  if (!update.all.pubmed & file.exists(paste(path, 'r/url_analysis2pubmed.rds', sep='/'))) {
    pmid<-readRDS(file=paste(path, 'r/analysis2pubmed.rds', sep='/')); 
    url1<-url1[!(names(url1) %in% names(pmid))];
  } else {
    pmid<-list();
  }
  if (length(url1) > 0) {
    pmid0<-list();
    for (i in 1:length(url1)) {
      print(i);
      pmid0[[i]]<-getPMID(url1[i]);
    }
    names(pmid0)<-names(url1);
    pmid<-c(pmid, pmid0);
  }
  pmid[sapply(pmid, length)==0]<-'-1';
  pmid<-pmid[rownames(ana)];
  
  ###############################################
  id2ana<-lapply(names(pmid), function(nm) {
    if (is.null(id2ana[[nm]][['PubMed']])) {
      append(id2ana[[nm]], list(PubMed=pmid[[nm]]), after=which(names(id2ana[[nm]])=='URL'));
    } else id2ana[[nm]];
  });
  names(id2ana)<-sapply(id2ana, function(x) x[[1]]);
  saveRDS(id2ana,file=paste(path, 'r/analysis_by_id.rds', sep='/'));
  
  ##################################################################
  ##################################################################
  
  
  #########
  # Study #
  #########
  
  ###############################################
  url0<-as.vector(std$URL);
  names(url0)<-rownames(std);
  
  # Get URL to CGI script that reports pubmed linked to each study
  if (!update.all.pubmed & file.exists(paste(path, 'r/url_study2pubmed.rds', sep='/'))) {
    url1<-readRDS(file=paste(path, 'r/url_study2pubmed.rds', sep='/'));
    url0<-url0[!(names(url0) %in% names(url1))];
  } else {
    url1<-c();
  }
  if (length(url0) > 0) {
    url.new<-sapply(url0, function(url) {
      html<-strsplit(getURL(url), '\n')[[1]];
      ln<-html[grep('initializeReferences', html)][1];
      ln<-strsplit(ln, '[\";]')[[1]];
      ln<-ln[grep('^initializeReferences', ln)][1]
      ln<-strsplit(ln, '[(\',)]')[[1]]; 
      print(ln[3]);   
      if (length(ln) !=7) '' else {
        paste("http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetReference.cgi?study_id=", ln[3], '&study_key=', ln[5], '&page_number=1', sep='');
      }
    });
    names(url.new)<-names(url0);
    url1<-c(url1, url.new);
  }
  url1<-url1[rownames(std)];
  saveRDS(url1, file=paste(path, 'r/url_study2pubmed.rds', sep='/'));
  
  ###############################################
  # Get Pubmed mapped to the study from a cgi script of dbGap
  url1<-readRDS(file=paste(path, 'r/url_study2pubmed.rds', sep='/'));
  if (!update.all.pubmed & file.exists(paste(path, 'r/url_study2pubmed.rds', sep='/'))) {
    pmid<-readRDS(file=paste(path, 'r/study2pubmed.rds', sep='/')); 
    url1<-url1[!(names(url1) %in% names(pmid))];
  } else {
    pmid<-list();
  }
  if (length(url1) > 0) {
    pmid0<-list();
    for (i in 1:length(url1)) {
      print(i);
      pmid0[[i]]<-getPMID(url1[i]);
    }
    names(pmid0)<-names(url1);
    pmid<-c(pmid, pmid0);
  }
  pmid[sapply(pmid, length)==0]<-'-1';
  pmid<-pmid[rownames(std)];
  
  ###############################################
  names(id2std)<-sapply(id2std, function(x) x$ID);
  id2std<-lapply(names(pmid), function(nm) {
    if (is.null(id2std[[nm]][['PubMed']])) {
      append(id2std[[nm]], list(PubMed=pmid[[nm]]), after=which(names(id2std[[nm]])=='URL'));
    } else id2std[[nm]];
  });
  names(id2std)<-sapply(id2std, function(x) x$ID);
  id2std<-id2std[rownames(std)];
  saveRDS(id2std,file=paste(path, 'r/study_by_id.rds', sep='/'));
  
  #################################################################
  # study/analysis to pubmed id as list
  id2ana<-readRDS(file=paste(path, 'r/analysis_by_id.rds', sep='/'));
  id2std<-readRDS(file=paste(path, 'r/study_by_id.rds', sep='/'));
  names(id2ana)<-sapply(id2ana, function(x) x$ID);
  names(id2std)<-sapply(id2std, function(x) x$ID);
  
  ana2pm<-lapply(id2ana, function(x) x$PubMed);
  std2pm<-lapply(id2std, function(x) x$PubMed);
  saveRDS(ana2pm,file=paste(path, 'r/analysis2pubmed.rds', sep='/'));
  saveRDS(std2pm,file=paste(path, 'r/study2pubmed.rds', sep='/'));
  
  #################################################################
  #################################################################
  pmid<-c(lapply(id2ana, function(x) x$PubMed), lapply(id2std, function(x) x$PubMed)); # All PubMed IDs
  pmid<-sort(unique(unlist(pmid)));
  pmid<-pmid[pmid!='-1' & pmid!='' & !is.na(pmid)];
  
  if (!update.all.pubmed & file.exists(paste(path, 'r/pubmed_downloaded.rds', sep='/'))) {
    pm<-readRDS(file=paste(path, 'r/pubmed_downloaded.rds', sep='/'));
    pmid0<-pmid[!(pmid %in% names(pm))];
    if (length(pmid0) > 0) {
      pm0<-GetPubMedAbstract(pmid0);
      #names(pm0)<-pmid0;
      pm<-append(pm, pm0);
    }
  } else {
    pm<-GetPubMedAbstract(pmid);
  }
  
  saveRDS(pm, file=paste(path, 'r/pubmed_downloaded.rds', sep='/'));
  pubmed<-GetPubMedFields(pm);
  saveRDS(pubmed, file=paste(path, 'r/pubmed.rds', sep='/'));
  
  ###################################################################################################################
  ###################################################################################################################
}