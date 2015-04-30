#!/usr/bin/env Rscript

# Copy all R function files in the dev/ folder to the R/ folder
fn<-dir('./dev', recursive=TRUE);
fn<-fn[grep('.r', fn, ignore.case=TRUE)];
fn0<-paste('./dev/', fn, sep=''); # source file names
fn<-sapply(strsplit(fn, '/'), function(x) x[length(x)]);
fn1<-paste('./R/', fn, sep='');
mt0<-file.info(fn0)[, 'mtime']; # compare modify time
mt1<-file.info(fn1)[, 'mtime'];
names(fn0)<-fn1;
fn0<-fn0[is.na(mt1) | mt0>mt1];
if (length(fn0)) {
	cat("copied", length(fn0), "function files to R/\n");
	file.copy(fn0, names(fn0), overwrite=TRUE);
}

# Export functions
fn<-paste('./R/', dir('./R/'), sep='');
fn<-fn[grep('.r$', fn, ignore.case=TRUE)];
if (length(fn)>0) {
	fnc<-lapply(fn, function(fn) {
		ln<-scan(fn, what='', flush=TRUE, sep='\n');
    ln<-ln[!grepl('^ ', ln)];
		ln<-gsub(' ', '', ln);
		ln<-ln[grep("<-function\\(", ln)];
		sapply(strsplit(ln, "<-function\\("), function(x) x[1]);
	});
	fnc<-sort(as.vector(unlist(fnc)));
  ln<-paste('export("', fnc, '");', sep='');
  writeLines(c('', ln, ''), file('./NAMESPACE'));
}
