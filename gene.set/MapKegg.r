##################################
# Mapping between KEGG databases #
##################################

######################################################################################################
######################################################################################################
# Use the KEGG dbget engine to map a ID to IDs of specified KEGG database
# Download and parse from an HTML page to get the mapping results
# The URL of the HTML page will be built as:
#     http://www.genome.jp/dbget-bin/get_linkdb?-t+<DB>+<Abbr>:<ID>
#         - <DB> is the name of the database to be mapped to
#         - <Abbr> is the abbrevation of the database defining the query ID
#         - <ID> is the query ID
# For example, "http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:hsa01100" maps pathway hsa01100 to genes
# Refer to this page for available databases, database names, name abbreviations: http://www.genome.jp/dbget/
MapKeggOne2One<-function(id, type, to, verbose=FALSE) {
  # id    A unique KEGG ID, gene, compound, drug, etc.
  # type  Type of the ID, use database abbreviation such as 'gn', 'path', and 'cpd', see http://www.genome.jp/dbget
  # to    The database that the id will be mapped to, use database name, such as 'genes', 'pathway', and 'compound'
  
  library(devtools);
  source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/extra/htmlUtilities.r"); # functions that remove HTML tags from a string
  
  # URL prefix
  url<-paste("http://www.genome.jp/dbget-bin/get_linkdb?-t+", to, '+', type, ':', id, sep='');
  
  download.file(url, 'tmp.html');
  html<-scan('tmp.html', sep='\n', what='', flush=TRUE);
  file.remove('tmp.html');
  
  if (verbose) {
    cat(url, '\n');
    cat(length(html), 'lines\n');
  }

  # Parsing HTML file
  ln<-html[(grep('^---', html)+1):(grep('integrated database retrieval system', html)-1)];
  ln<-gsub('</a> ', '\t</a>', ln);
  ln<-CleanHtmlTags(ln);
  ln<-gsub('\t +', '\t', ln);
  ln<-gsub('\t$', '\t ', ln);
  
  # Turn to named vectors
  x<-do.call('rbind', strsplit(ln, '\t'));
  mp<-x[,2];
  names(mp)<-x[,1];
  mp[mp==' ' | is.na(mp)]<-'';
  
  mp;
}
  
######################################################################################################
######################################################################################################
# Use the KEGG dbget engine to map a ID to IDs of all KEGG database
# Download and parse from an HTML page to get the mapping results
# The URL of the HTML page will be built as:
#     http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+<Abbr>:<ID> 
#         - <Abbr> is the abbrevation of the database defining the query ID
#         - <ID> is the query ID
# For example, "http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+path:hsa01100" maps pathway hsa01100 to IDs of all databases
# Refer to this page for available databases and name abbreviations: http://www.genome.jp/dbget/
MapKeggOne2All<-function(id, type, verbose=FALSE) {
  # id    A unique KEGG ID; gene, compound, drug, etc.
  # type  Type of the ID; such as 'gn', 'cpd', etc. see http://www.genome.jp/dbget
  
  #library(RCurl);
  library(devtools);
  source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/extra/htmlUtilities.r"); # functions that remove HTML tags from a string
  
  # URL prefix
  prefix<-"http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+";
  
  types<-c(
    'pathway' = 'path', 'path' = 'path',
    'brite' = 'br', 'br' = 'br',
    'module' = 'md', 'md' = 'md',
    'enzyme' = 'ec', 'ec' = 'ec',
    'reaction' = 'rn', 'rn' = 'rn',
    'rpair' = 'rp', 'rp' = 'rp', 
    'rclass' = 'rc', 'rc' = 'rc', 
    'orthology' = 'ko', 'ko' = 'ko', 
    'gene' = 'gn', 'gn' = 'gn', 
    'disease' = 'ds', 'ds' = 'ds',
    'compound' = 'cpd', 'cpd' = 'cpd',
    'glycan' = 'gl', 'gl' = 'gl',
    'drug' = 'dr', 'dr' = 'dr',
    'dgroup' ='dg', 'dg' = 'dg',
    'environ' = 'ev', 'ev' = 'ev'
  )
  tp<-types[tolower(type)];
  if (is.na(tp)) tp<-type;
  
  if (verbose) cat(tp, '=', id, '\n');
  
  # Read in HTML file from URL
  url<-paste(prefix, tp, ':', id, sep='');
  #html<-strsplit(getURL(url), '\n')[[1]];
  download.file(url, 'tmp.html');
  html<-scan('tmp.html', sep='\n', what='', flush=TRUE);
  file.remove('tmp.html');
  
  if (verbose) {
    cat(url, '\n');
    cat(length(html), 'lines\n');
  }
  ind.db<-grep('^<span ', html); # lines specifying source databases 
  
  if (length(ind.db) == 0) list() else {
    html<-gsub('</a> ', '\t</a>', html);
    ln<-CleanHtmlTags(html); # Remove HTML tags
    ln<-gsub('\t +', '\t', ln);
    ln<-gsub('\t$', '\t ', ln);
    
    # Mapping
    nm.db<-ln[ind.db];
    n<-c(ind.db[-1]-ind.db[-length(ind.db)], length(ln)-ind.db[length(ind.db)]+1);
    nm<-rep(c('', nm.db), c(ind.db[1]-1, n));
    mp<-split(ln, nm)[-1];
    mp<-lapply(mp, function(mp) mp[-1]);
    mp<-lapply(mp, function(mp) mp[mp!='' & !is.na(mp)]);
    mp<-mp[sapply(mp, length)>0];
    
    # Turn to named vectors
    mp<-lapply(mp, function(mp) {
      x<-do.call('rbind', strsplit(mp, '\t'));
      y<-x[,2];
      names(y)<-x[,1];
      y[y==' ' | is.na(y)]<-'';
      y;
    })
  
    mp;
  }
}