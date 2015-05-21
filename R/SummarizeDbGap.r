################################################################################
# Summarize metadata infomration of dbGaP studies
SummarizeDbGap<-function(meta, path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'), update.all.study=FALSE) {
  # meta              Analysis metadata table prepared by the DownloadDbGap function
  # path              Path to output folder
  # update.all.study  If TRUE, re-download information of all studies from source; Just download the new ones if FASLE
  
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
      #print(f);
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
    st<-lapply(st.id, function(x) GetDbGapStudy(x)); # get study info from dbGAP web page
    names(st)<-st.id;
    st<-st[order(names(st))];
    saveRDS(st, file=fn.dbgap);
  }
  
  ###############################################
  # analysis anotation table
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
