########## Functions retrieving data from PubMed ########## 

MapPMID2Title <- function(pmid, set.size=200) {
  getSet <- function(pmid) {
    url <- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=';
    # url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=json&rettype=title&id=";
    url <- paste(url, paste(pmid, collapse=','), sep='');
    jsn <- jsonlite::fromJSON(url)$result;
    # lns <- readLines(url);
    # lns <- paste(lns, collapse='\n'); 
    # lns <- strsplit(lns, '\n\n')[[1]];
    # ttl <- lns[seq(2, length(lns), 4)]; 
    pid <- jsn$uids;
    ttl <- sapply(jsn[pid], function(j) j$title); 
    names(ttl) <- pid;
    ttl;
  }

  pmid<-pmid[!is.na(pmid) & pmid!=''];
  pmid<-sort(unique(pmid));
  if (length(pmid) > 0) {
    n<-length(pmid);
    n.set<-ceiling(n/set.size);
    sets<-lapply(1:n.set, function(i) {
      cat(set.size*i, '\n');
      getSet(pmid[(i*set.size-set.size+1):min(i*set.size, length(pmid))]);
    });
    do.call('c', sets); # return list
  } else {
    cat("No valid PMID provided\n");
    list();
  }
}

# Retrieve PubMed abstract and other information using the efetch util of NCBI and parse the XML document into a list
GetPubMedAbstract<-function(pmid, set.size=200) {
  # pmid      One or a set of PMID ids
  # set.size  Maximum number of PMID articles to be retrieved together at one time. Set size over 500 might cause URL length over limit
  
  library(XML);
  
  # Retrieve up to 500 pubmed IDs
  getSet<-function(pmid) {
    cat('Retrieve', length(pmid), 'PubMed articles\n');
    # build URL
    url<-"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&rettype=abstract&id=";
    url<-paste(url, paste(pmid, collapse=','), sep='');
    
    # retrieve xml document
    xml.doc<-xmlTreeParse(url, useInternalNodes=TRUE);
    xml.top<-xmlRoot(xml.doc);
    
    # convert to list
    xml2list<-xmlToList(xml.top);
    names(xml2list)<-sapply(xml2list, function(x) x$MedlineCitation$PMID$text);
    
    xml2list;
  }
  
  pmid<-pmid[!is.na(pmid) & pmid!=''];
  pmid<-sort(unique(pmid));
  if (length(pmid) > 0) {
    n<-length(pmid);
    n.set<-ceiling(n/set.size);
    sets<-lapply(1:n.set, function(i) getSet(pmid[(i*set.size-set.size+1):min(i*set.size, length(pmid))]));
    do.call('c', sets); # return list
  } else {
    cat("No valid PMID provided\n");
    list();
  }
}

# Get fields of Pubmed entries from output from the GetPubMedAbstract function, and return a data.frame
GetPubMedFields<-function(articles, fields=c('Title', 'Journal', 'Date', 'URL', 'Abstract')) {
  # articles    A list of Pubmed articles, each including nested fields similar to nodes of XML document. Usually generated by the GetPubMedAbstract function
  # fields      The name of the fields to be retrieved and column names of the output data.frame
  
  fields<-unique(c('Title', fields)); # must include article title
  
  # field name mapping
  mp<-c(
    'title' = "MedlineCitation.Article.ArticleTitle",
    'journal' = "MedlineCitation.Article.Journal.Title",
    'issn' = "MedlineCitation.Article.Journal.ISSN.text",
    'volumn' = "MedlineCitation.Article.Journal.JournalIssue.Volume",                            
    'issue' = "MedlineCitation.Article.Journal.JournalIssue.Issue",
    'year' = "MedlineCitation.DateCreated.Year",
    'month' = "MedlineCitation.DateCreated.Month",
    'day' = "MedlineCitation.DateCreated.Day",
    'lastname' = "MedlineCitation.Article.AuthorList.Author.LastName",
    'initials' = "MedlineCitation.Article.AuthorList.Author.Initials",
    'affiliation' = "MedlineCitation.Article.AuthorList.Author.AffiliationInfo.Affiliation"
  );
  # Field names of abstract
  abs.nm = c("MedlineCitation.Article.Abstract.AbstractText", "MedlineCitation.Article.Abstract.AbstractText.text");
  
  # get article PMID and corresponding URL
  art<-lapply(articles, unlist);
  pmid<-sapply(art, function(art) art['MedlineCitation.PMID.text']);
  url<-paste('http://www.ncbi.nlm.nih.gov/pubmed/?term=', pmid, sep='');
  
  # Create the output table
  tbl<-t(sapply(art, function(art) art[mp]));
  colnames(tbl)<-names(mp);
  rownames(tbl)<-pmid;
  tbl<-cbind(tbl, url=url, date=paste(tbl[, 'year'], tbl[, 'month'], tbl[, 'day'], sep=''), 
             author=paste(tbl[, 'lastname'], ', ', tbl[, 'initials'], sep=''));
  tbl<-cbind(tbl, abstract=sapply(art, function(art) paste(art[names(art) %in% abs.nm], collapse=' ')));
  tbl<-tbl[!duplicated(rownames(tbl)), , drop=FALSE];
  
  # select column in output
  colnames(tbl)<-tolower(colnames(tbl));
  fld<-tolower(fields);
  names(fld)<-fields;
  fld<-fld[fld %in% colnames(tbl)];
  tbl<-tbl[, fld, drop=FALSE];
  colnames(tbl)<-names(fld);

  data.frame(tbl, stringsAsFactors=FALSE);
}