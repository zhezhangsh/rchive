# Use the GEOquery::getGEO function to retrieve a GEO data series and parse the essential data of the series

ProcessGEO<-function(gse.id, path.out, affy.supp=TRUE, gsm.ids=c()) {  
  # gsm.ids     A subset samples of the same GSE, must from the same platform
  
  # download GEO data sets
  if (length(gsm.ids)>0) path<-ParseGSM(gsm.ids, paste(path.out, gse.id, sep='/'), affy.supp) else 
    path<-ParseGSE(gse.id, path.out, affy.supp);
  
  if (affy.supp) {
    path.smp<-paste(path, gse.id, sep='/');
    if (length(gsm.ids) == 0) untar(paste(path.smp, '/', gse.id, '_RAW.tar', sep=''), exdir=path.smp);

    fn.cel<-dir(path.smp);
    fn.cel<-paste(path.smp, fn.cel[grep('.CEL.gz$', fn.cel, ignore.case=TRUE)], sep='/');
    
    raw<-LoadAffyCel(fn.cel);
    raw@cdfName<-InstallBrainarray(raw@cdfName);
    if (!identical(NA, raw@cdfName)) {
      expr<-exprs(rma(raw));
      expr<-expr[grep('_at$', rownames(expr)), , drop=FALSE];
      rownames(expr)<-sub('_at$', '', rownames(expr));
      cnm<-sub('.CEL.gz$', '', sampleNames(raw));
      cnm<-sapply(strsplit(cnm, '/'), function(x) x[length(x)]);
      colnames(expr)<-cnm;
      saveRDS(expr, file=paste(path.out, '/expr_', gse.id, '.rds', sep=''));      
     }
  }
  
  smp<-eval(parse(text=load(paste(path, 'phen.rdata', sep='/'))));
  #smp<-smp[colnames(expr), , drop=FALSE];
  smp<-cbind(smp[, c('title', 'taxid_ch1')], 
             smp[, c(grep('characteristics_ch1', colnames(smp)), grep('^description', colnames(smp)))]);
  saveRDS(smp, file=paste(path.out, '/sample_', gse.id, '.rds', sep=''));      
  
  path;
}

# Load a set of CEL files
LoadAffyCel<-function(fn) {
  fn<-fn[file.exists(fn)];
  if (length(fn) == 0) NA else {
    heads<-sapply(fn, function(f) affyio::read.celfile.header(f)[[1]]);
    if (length(unique(heads)) > 1) NA else {
      hd<-affyio::read.celfile.header(fn[1]);
      cdfname<-hd[[1]];
      dim.int<-hd[[2]];
      exprs<-affyio::read_abatch(fn, FALSE, FALSE, FALSE, cdfname, dim.int[c(1, 2)], TRUE);
      raw<-new("AffyBatch", exprs = exprs, cdfName = cdfname, nrow = dim.int[2], ncol = dim.int[1],
               annotation = cleancdfname(cdfname, addcdf = FALSE));

    }
  }
}

# Install a BrainArray package;
InstallBrainarray<-function(affy.name, type='entrezg', version='19.0.0') {
  url="http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF"
  u<-paste(url, '/', version, '/', type, '.asp', sep='');
  if(!RCurl::url.exists(u)) NA else {
    tbls<-XML::readHTMLTable(u);
    tbl<-tbls[[length(tbls)]];
    
    all<-tolower(tbl[, 3]);
    nm.old<-cleancdfname(affy.name, addcdf = FALSE)
    ind<-which(all == nm.old);
    if (length(ind) == 0) NA else {
      nm.new<-tbl[ind, 5];
      nm.new<-tolower(gsub('[-_]', '', nm.new));
      urls<-paste("http://mbni.org/customcdf/", version, '/', type, '.download/', nm.new, 
                  c('cdf', 'probe', '.db'), '_', version, '.tar.gz', sep='');
      urls<-urls[sapply(urls, RCurl::url.exists)];
      if (length(urls) == 0) NA else {
        fn<-sapply(urls, install_url);
        nm.new;
      }
    }
  } 
}

ParseGSM<-function(ids, destdir=getwd(), getSupp=TRUE) {
  library(GEOquery);
  
  cat('Parsing', length(ids), 'samples\n');
  
  if (!file.exists(destdir)) dir.create(destdir);
  gsm<-lapply(ids, function(id) getGEO(id, destdir=destdir, AnnotGPL=FALSE));
  
  phen<-t(sapply(gsm, function(sm) unlist(sm@header)));
  phen<-data.frame(phen, stringsAsFactors=FALSE);
  rownames(phen)<-ids;
  
  plt.id<-gsm[[1]]@header$platform_id;
  anno<-getGEO(plt.id);
  anno<-anno@dataTable@table;
  #rownames(anno)<-anno[, 'ID'];
  #anno<-anno[, -1];
  
  expr<-sapply(gsm, function(sm) sm@dataTable@table);
  
  save(expr, file=paste(destdir, '/expr.rdata', sep=''));
  save(anno, file=paste(destdir, '/anno.rdata', sep=''));
  save(phen, file=paste(destdir, '/phen.rdata', sep=''));
  
  #write.table(cbind(ID=rownames(expr), expr), file=paste(f, '/expr.txt', sep=''), sep='\t', qu=FALSE, row=FALSE); 
  write.table(anno, file=paste(destdir, '/anno.txt', sep=''), sep='\t', qu=FALSE, row=FALSE);
  write.table(phen, file=paste(destdir, '/phen.txt', sep=''), sep='\t', qu=FALSE, row=FALSE); 
  
  if (getSupp) {
    path.smp<-paste(destdir, gsm[[1]]@header$series_id, sep='/'); print(path.smp)
    if (!file.exists(path.smp)) dir.create(path.smp);
    supp<-sapply(ids, function(id) getGEOSuppFiles(id, TRUE, path.smp));
    fn<-dir(path.smp, rec=TRUE);
    file.rename(paste(path.smp, fn, sep='/'), paste(path.smp, sapply(strsplit(fn, '/'), function(x) x[length(x)]), sep='/'));
  }
  
  destdir;
}

ParseGSE<-function(id, destdir=getwd(), getSupp=TRUE) { 
  # id        GEO data series ID
  # dir       Location where a subdirectory will be created to store the data
  # getSupp  Download supplemental files if TRUE
  
  library(GEOquery);
  
  cat('Parsing GEO series', id, '\n');
  
  if (!file.exists(destdir)) dir.create(destdir);
  gse<-getGEO(id[1], destdir=destdir, AnnotGPL=FALSE)[[1]];
  
  gpl<-gse@annotation; 
  f<-paste(destdir, '/', id, '-', gpl, sep='');
  if(!file.exists(f)) dir.create(f);
  
  expr<-exprs(gse);
  anno<-featureData(gse)@data;
  phen<-phenoData(gse)@data
  
  save(expr, file=paste(f, '/expr.rdata', sep=''));
  save(anno, file=paste(f, '/anno.rdata', sep=''));
  save(phen, file=paste(f, '/phen.rdata', sep=''));
  
  write.table(cbind(ID=rownames(expr), expr), file=paste(f, '/expr.txt', sep=''), sep='\t', qu=FALSE, row=FALSE); 
  write.table(anno, file=paste(f, '/anno.txt', sep=''), sep='\t', qu=FALSE, row=FALSE);
  write.table(phen, file=paste(f, '/phen.txt', sep=''), sep='\t', qu=FALSE, row=FALSE); 
  
  if (getSupp) {
    getGEOSuppFiles(id, TRUE, f);
  }
  
  f;
}