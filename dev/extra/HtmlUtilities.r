### Utility functions used to process HTML files

########################################################################################################################
## clean up HTML tags in a line
#' @export
CleanHtmlTags<-function(l) {
  l<-gsub("<.*?>", " ", l); # remove open tags
  l<-gsub("</.*?>", " ", l); # remove close tags
  l<-gsub(" +", " ", l); # replace multiple spaces with one
  l<-sub('^ ', '', l); # remove space at the begining of the line
  l<-sub(' $', '', l); # remove space at the end of the line
  l;
}
########################################################################################################################librar
