### Utility functions used to process HTML files

########################################################################################################################
## clean up HTML tags in a line
#' @export
CleanHtmlTags<-function(l, remove.empty=TRUE) {
  l<-gsub("<.*?>", " ", l); # remove open tags
  l<-gsub("</.*?>", " ", l); # remove close tags
  l<-gsub(" +", " ", l); # replace multiple spaces with one
  l<-sub('^ ', '', l); # remove space at the begining of the line
  l<-sub(' $', '', l); # remove space at the end of the line
  
  if (remove.empty) l<-l[l!='' & !is.na(l)];
  l;
}

########################################################################################################################librar
# Download a set of URLs in parellel
GetSetURL<-function(url, max.try=10, split.lines=TRUE, clean.tags=TRUE, select.lines=NA) {
  # url             One or a set of URLs
  # max.try         Max number of tries if no more new pages were downloaded
  # split.lines     Split lines on each page if TRUE, or <select.lines> is a non-empty string
  # clean.tags      Clean up all HTML tags if TRUE
  # selected.lines  Select the lines having the given string, using regular expression
  
  library(RCurl);
  
  # empty pages
  pgs<-rep('', length(url));
  
  # number of tries
  n.try<-0;
  
  # number of un-downloaded pages
  n.pgs<-length(pgs[pgs=='']);
  
  #################################################
  # Get 100 pages a time
  get.by.100<-function(url) {
    url<-split(url, rep(1:ceiling(length(url)/100), each=100, length.out=length(url)));
    pgs<-lapply(url, getURL);
    unlist(pgs, use.names=FALSE);
  }
  #################################################
  
  while(length(pgs[is.na(pgs) | pgs==''])>0 & n.try<=max.try) { 
    # download un-downloaded pages
    pgs[is.na(pgs) | pgs=='']<-get.by.100(url[is.na(pgs) | pgs=='']);
    
    if (n.pgs == length(pgs[is.na(pgs) | pgs==''])) { # no new downloads achieved
      n.try<-n.try+1;
    } else {
      n.try<-0;
      n.pgs<-length(pgs[is.na(pgs) | pgs=='']);
      if (n.pgs > 0) cat(n.pgs, 'pages to be downloaded\n');
    }
  }
  
  if (length(pgs[is.na(pgs) | pgs=='']) > 0) warning(n.pgs, ' pages have no contents');
  
  # Clean up HTML tags
  #if (clean.tags) pgs<-lapply(pgs, CleanHtmlTags);
  
  # Select lines including given sub-string if specified
  select.lines<-select.lines[select.lines!='' & !is.na(select.lines)]
  if (split.lines | length(select.lines)>0) {
    # split each page into lines
    pgs<-lapply(pgs, function(pg) strsplit(pg, '\n')[[1]]);
    if (clean.tags) pgs<-lapply(pgs, CleanHtmlTags); # Clean up HTML tags
    if (length(select.lines) == 1) {
      pgs<-lapply(pgs, function(pg) pg[grep(select.lines[1], pg)]);
    } else if (length(select.lines) > 1) {
      pgs<-lapply(pgs, function(pg) lapply(select.lines, function(s) pg[grep(s, pg)]));
    }
  } else if (clean.tags) pgs<-lapply(pgs, CleanHtmlTags);
  
  names(pgs)<-url;
  
  pgs;
}