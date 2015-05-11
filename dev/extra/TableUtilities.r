### Utility functions used to deal with matrix or data.frame

########################################################################################################################
# Load a specified column from a set of files, and create a joint table with each column corresponds to a 
LoadColumn<-function(fns, cnm, row=0, file.type=c('rds', 'rdata'), default.value=NA, collapse=c('none', 'first', 'last', 'mean', 'max', 'min'), is.numeric=FALSE, as.data.frame=FALSE) {
  # fns           A set of input file names; each file must be a R table saved as .rds or .rdata
  # cnm           Column name(s); use the first if more than one matches
  # row           The row name of the output table; the original rownames by default; Otherwise, specify the column name or index of the original table
  # file.type     Type of input file
  # default.value Default of initial and missing values
  # collapse      Options to collapse values having the same name; the options 'mean', 'max', 'min' only apply when <is.numeric> is TRUE
  # is.numeric    Make a numeric table if TRUE
  # as.data.frame If TRUE format output as a data.frame; as a matrix otherwise
  
  fns<-fns[file.exists(fns)];
  
  if (length(fns) == 0) out<-matrix(nr=0, nc=0) else {
    #################################################
    columns<-lapply(fns, function(fn) {
      # load file based on file type
      if (file.type[1]=='rds') tbl<-readRDS(fn) else tbl<-eval(parse(text=load(fn)));
      
      # Get columns
      c<-tbl[, tolower(colnames(tbl)) %in% tolower(cnm), drop=FALSE];
      if (is.numeric(row)) { # column index
        if (row[1]<=0 | row[1]>ncol(tbl)) r<-matrix(rownames(tbl), nc=1) else r<-tbl[, round(row[1]), drop=FALSE];
      } else { # column name
        r<-tbl[, tolower(colnames(tbl)) %in% tolower(row), drop=FALSE];
      }
      
      if (ncol(c)==0 | ncol(r)==0) c() else {  # Return empty vector if no matching columns
        if (is.numeric) c<-as.numeric(as.vector(c[,1])) else c<-as.vector(c[,1]);
        r<-as.vector(r[,1]); 
        names(c)<-r;
        c<-c[!is.na(c) & !is.na(r)];
        # collapse redundancy
        if (length(c) == 0) c() else {
          mp<-split(as.vector(c), names(c));
          mp<-lapply(mp, unique);
          if (collapse[1]=='first') {
            v<-sapply(mp, function(mp) mp[1]); 
          } else if (collapse[1]=='last') {
            v<-sapply(mp, function(mp) mp[length(mp)]); 
          } else if (collapse[1]=='mean' & is.numeric) {
            v<-sapply(mp, mean); 
          } else if (collapse[1]=='max' & is.numeric) {
            v<-sapply(mp, max); 
          } else if (collapse[1]=='min' & is.numeric) {
            v<-sapply(mp, min); 
          } else {
            v<-unlist(mp, use.names=FALSE);
            names(v)<-rep(names(mp), sapply(mp, length));
          }
          v;
        }
      }      
    }) # end of columns<-
    #################################################
    
    # making the table
    n<-sapply(columns, length);
    if (max(n) == 0) out<-matrix(nr=0, nc=0) else {
      row.id<-lapply(columns, names);
      row.id<-sort(unique(unlist(row.id, use.names=FALSE)));
      col.id<-names(columns);
      if (is.null(col.id)) col.id<-1:length(columns);
      col.id[is.na(col.id)]<-(1:length(fns))[is.na(col.id)];
      out<-matrix(default.value, nr=length(row.id), nc=length(col.id), dimnames=list(row.id, col.id));
      for (i in 1:length(columns)) out[names(columns[[i]]), i]<-as.vector(columns[[i]]);
    }
    out;
  }
  
  if (as.data.frame) data.frame(out, stringsAsFactors=FALSE) else out;
}


########################################################################################################################
# Merge a set of table to create a joint one with union sets of row and column names
MergeTable<-function(tbls, default.value=NA, use.first=FALSE, as.data.frame=FALSE) {
  # tbls            A list of tables (matrix or data.frame) to be merged
  # default.value Default of initial and missing values
  # use.first       When there is redundance (same row and column in multiple tables), use the first occurance if TRUE; use the last otherwise
  # as.data.frame   Whether the merged table should be matrix or data.frame
  
  # unique, union set of row names
  rnm<-sort(unique(unlist(lapply(tbls, rownames), use.names=FALSE)));

  # unique, union set of column names
  cnm<-sort(unique(unlist(lapply(tbls, colnames), use.names=FALSE)));
  
  out<-matrix(default.value, nr=length(rnm), nc=length(cnm), dimnames=list(rnm, cnm));
  
  if (use.first) tbls<-rev(tbls);
  
  for (i in 1:length(tbls)) out[rownames(tbls[[i]]), colnames(tbls[[i]])]<-tbls[[i]];
  
  if (as.data.frame) as.data.frame(out, stringsAsFactors=FALSE) else out;  
}
