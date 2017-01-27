yml<-yaml::yaml.load_file(paste(Sys.getenv("RCHIVE_HOME"), "source/script/update/variant/LoadAnnoVar.yml", sep='/')); 

path<-paste(Sys.getenv("RCHIVE_HOME"), 'data/variant/public/annovar/r', sep='/'); 
if (!file.exists(path)) dir.create(path); 

library(GenomicRanges)
library(rtracklayer)

fn<-lapply(names(yml$file), function(fn) {
  cat(fn, '\n');
  
  f<-yml$file[[fn]]; 
  tbl<-read.table(f, sep='\t', header = FALSE, comment.char = '#', stringsAsFactors = FALSE);
  hd<-scan(f, what='', flush=TRUE, sep='\n', nlines = 1); 
  if (grepl('^#', hd)) {
    sub('^#', '', hd); 
    cnm<-strsplit(hd, '\t')[[1]];
    colnames(tbl)<-cnm; 
  }
  if (class(tbl[, 3])!='integer') tbl<-cbind(tbl[, 1:2], tbl[, 2], tbl[, 3:ncol(tbl)]); 
  colnames(tbl)[1:3]<-c('seqnames', 'start', 'end');
  gr<-as(tbl, 'GRanges'); 
  names(gr)<-1:length(gr) ;
  saveRDS(gr, paste(path, paste(fn, '.rds', sep=''), sep='/')); 
  fn;
}); 


##################################################################################################
# Save run
tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/variant/LoadAnnoVar.r', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/variant/log/', tm, '_LoadAnnoVar.r' , sep='');
file.copy(fn0, fn1, overwrite = TRUE); 
fn0<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/variant/LoadAnnoVar.yml', sep='');
fn1<-paste(Sys.getenv("RCHIVE_HOME"), '/source/script/update/variant/log/', tm, '_LoadAnnoVar.yml' , sep='');
file.copy(fn0, fn1, overwrite = TRUE); 
