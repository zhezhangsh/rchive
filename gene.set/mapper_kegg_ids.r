# Given a KEGG ID, map it to entities of other KEGG databases, by directly downloading and parsing an HTML file
mapper_kegg_ids<-function(id, type) {
  # id    A unique KEGG ID; gene, compound, drug, etc.
  # type  Type of the ID; such as 'gene', 'compound', etc.
  
  library(RCurl);
    
  # URL prefix
  prefix<-"http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+";
    
  mp<-c(
    'gene' = 'gn', 'gn' = 'gn', 
    'compound' = 'cpd', 'cpd' = 'cpd',
    'drug' = 'dr', 'dr' = 'dr'
  )
  tp<-mp[tolower(type)];
  
  if (is.na(prefix)) list() else {
    url<-paste(prefix, tp, ':', id, sep='');
    html<-strsplit(getURL(url), '\n')[[1]];
  }
}