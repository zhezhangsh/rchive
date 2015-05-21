################################################################################
# Associate dbGaP analysis and study with Pubmed articles
AddDbGapPubMed<-function(path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gwas/public/dbgap', sep='/'), update.all.pubmed=FALSE) {
  # path              Path to the existing dbGaP data previously parsed 
  # update.all.pubmed If TRUE, update all pubmed-related information; otherwise, using previously saved information when it's available
  
  ###################################################################################################################
  ###################################################################################################################
  # Now, 
  std<-readRDS(file=paste(path, 'r/study.rds', sep='/'));
  ana<-readRDS(file=paste(path, 'r/analysis.rds', sep='/'));
  id2std<-readRDS(file=paste(path, 'r/study_by_id.rds', sep='/'));
  id2ana<-readRDS(file=paste(path, 'r/analysis_by_id.rds', sep='/'));
  
  library(RCurl);
  library(rchive);
  
  ##################################################################
  # get one page of Pubmed IDs based on given URL
  getPMID<-function(urls) {
    url<-unlist(urls, use.names=FALSE);
    id<-rep(names(urls), sapply(urls, length));
    
    pgs<-GetSetURL(url, clean.tags=FALSE, select.lines='db=pubmed'); # Download from dbGaP CGI page
    pmid<-lapply(pgs, function(ln) sapply(strsplit(ln, '[=\\[]'), function(ln) ln[length(ln)-1])); # parse page to get PubMed IDs
    
    # Group PubMed IDs
    mp<-split(pmid, id);
    mp<-lapply(mp, function(x) {
      ids<-unlist(x, use.names=FALSE);
      if (length(ids) == 0) '-1' else sort(unique(ids));
    });
    
    mp[names(urls)];
  }
  ##################################################################
  
  ########################
  # Analysis to pubmed
  ########################
  
  ###############################################
  ana<-readRDS(file=paste(path, 'r/analysis.rds', sep='/'));
  url0<-as.vector(ana$URL); # URL to each analysis
  names(url0)<-rownames(ana);
  
  # Get URL to CGI script that reports pubmed linked to each analysis
  if (!update.all.pubmed & file.exists(paste(path, 'r/url_analysis2pubmed.rds', sep='/'))) {
    url1<-readRDS(file=paste(path, 'r/url_analysis2pubmed.rds', sep='/'));
    url0<-url0[!(names(url0) %in% names(url1)) | is.na(url0) | url0==''];
  } else {
    url1<-c();
  }
  if (length(url0) > 0) {
    # Get URL of new pubmed 
    url.new<-paste("http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetAnalysisReference.cgi?pha=", 
                   as.numeric(sub('pha', '', names(url0))), '&version=1&page_number=1', sep=''); 
    names(url.new)<-names(url0);
    url1<-c(url1, url.new);
  }
  ###
  saveRDS(url1[rownames(ana)], file=paste(path, 'r/url_analysis2pubmed.rds', sep='/'));
  ###
  
  ###############################################
  # Get Pubmed mapped to the analysis from a cgi script of dbGap
  url1<-readRDS(file=paste(path, 'r/url_analysis2pubmed.rds', sep='/'));
  if (!update.all.pubmed & file.exists(paste(path, 'r/analysis2pubmed.rds', sep='/'))) {
    pmid<-readRDS(file=paste(path, 'r/analysis2pubmed.rds', sep='/')); 
    url1<-url1[!(names(url1) %in% names(pmid)) & !is.na(url1) & url1!=''];
  } else {
    pmid<-list();
  }
  if (length(url1) > 0) {
    pmid0<-getPMID(url1);
    pmid<-c(pmid, pmid0);
  }
  pmid<-pmid[rownames(ana)];
  
  ###############################################
  id2ana[1:length(id2ana)]<-lapply(names(id2ana), function(nm) {
    if (!is.null(pmid[[nm]])) {
      if (is.null(id2ana[[nm]][['PubMed']])) {
        append(id2ana[[nm]], list(PubMed=pmid[[nm]]), after=which(names(id2ana[[nm]])=='URL'));
      } else id2ana[[nm]][['PubMed']]<-pmid[[nm]];
    }
  });
  ###
  saveRDS(id2ana,file=paste(path, 'r/analysis_by_id.rds', sep='/'));
  saveRDS(lapply(id2ana, function(x) x$PubMed), file=paste(path, 'r/analysis2pubmed.rds', sep='/'));
  
  ###
  
  ##################################################################
  ##################################################################
  
  
  ########################
  # Study to pubmed
  ########################
  
  ###############################################
  std<-readRDS(file=paste(path, 'r/study.rds', sep='/'));
  url0<-as.vector(std$URL);
  names(url0)<-rownames(std);
  
  # Get URL to CGI script that reports pubmed linked to each study
  if (!update.all.pubmed & file.exists(paste(path, 'r/url_study2pubmed.rds', sep='/'))) {
    url1<-readRDS(file=paste(path, 'r/url_study2pubmed.rds', sep='/'));
    url0<-url0[!(names(url0) %in% names(url1)) | is.na(url0) | url0==''];
  } else {
    url1<-c();
  }
  if (length(url0) > 0) {
    pgs<-GetSetURL(url0, select.lines="initializeReferences", clean.tags=FALSE);
    lns<-sapply(strsplit(sapply(pgs, function(x) x[1]), '[\";]'), function(x) x[grep('^initializeReferences', x)][1]);
    lns<-strsplit(lns, '[(\',)]');
    url.new<-lapply(lns, function(ln) if (length(ln) !=7) list() else {
      url<-paste("http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetReference.cgi?study_id=", ln[3], '&study_key=', ln[5], sep='');
      ln<-GetSetURL(url, select.lines="There (is)|(are) [0-9]+ selected publications? related to this study.")[[1]];
      n<-as.numeric(strsplit(ln, ' ')[[1]][3]);
      paste(url, '&page_number=', 1:ceiling(n/10), sep='');
    });
    names(url.new)<-names(url0);
    url1<-c(url1, url.new);
  }
  ###
  saveRDS(url1[rownames(std)], file=paste(path, 'r/url_study2pubmed.rds', sep='/'));
  ###
  
  ###############################################
  # Get Pubmed mapped to the study from a cgi script of dbGap
  url1<-readRDS(file=paste(path, 'r/url_study2pubmed.rds', sep='/'));
  if (!update.all.pubmed & file.exists(paste(path, 'r/study2pubmed.rds', sep='/'))) {
    pmid<-readRDS(file=paste(path, 'r/study2pubmed.rds', sep='/')); 
    url1<-url1[!(names(url1) %in% names(pmid)) & !is.na(url1) & url1!=''];
  } else {
    pmid<-list();
  }
  # Retrieve new PubMed IDs from dbGaP
  if (length(url1) > 0) {
    pmid0<-getPMID(url1);
    pmid<-c(pmid, pmid0);
  }
  pmid<-pmid[rownames(std)];
  
  ###############################################
  names(id2std)<-sapply(id2std, function(x) x$ID);
  id2std[1:length(id2std)]<-lapply(names(id2std), function(nm) {
    if (!is.null(pmid[[nm]])) {
      if (is.null(id2std[[nm]][['PubMed']])) {
        append(id2std[[nm]], list(PubMed=pmid[[nm]]), after=which(names(id2std[[nm]])=='URL'));
      } else id2std[[nm]];
    }    
  });
  ###
  saveRDS(id2std,file=paste(path, 'r/study_by_id.rds', sep='/'));
  saveRDS(lapply(id2std, function(x) x$PubMed), file=paste(path, 'r/study2pubmed.rds', sep='/'));
  ###
  
  #################################################################
  #################################################################
  # Retrieve metadata of pubmed articles
  pmid<-c(lapply(id2ana, function(x) x$PubMed), lapply(id2std, function(x) x$PubMed)); # All PubMed IDs
  pmid<-sort(unique(unlist(pmid)));
  pmid<-pmid[pmid!='-1' & pmid!='' & !is.na(pmid)];
  
  if (!update.all.pubmed & file.exists(paste(path, 'r/pubmed_downloaded.rds', sep='/'))) {
    pm<-readRDS(file=paste(path, 'r/pubmed_downloaded.rds', sep='/'));
    pmid0<-pmid[!(pmid %in% names(pm))];
    if (length(pmid0) > 0) {
      pm0<-GetPubMedAbstract(pmid0);
      pm<-c(pm, pm0);
    }
  } else {
    pm<-GetPubMedAbstract(pmid);
  }
  
  saveRDS(pm, file=paste(path, 'r/pubmed_downloaded.rds', sep='/'));
  pubmed<-GetPubMedFields(pm);
  ###
  saveRDS(pubmed, file=paste(path, 'r/pubmed.rds', sep='/'));
  ###
}

