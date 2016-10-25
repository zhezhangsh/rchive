MergeEncodeReadCount <- function(fns) {
  require(RoCA); 
  
  fns <- fns[file.exists(fns)]; 
  
  d <- lapply(fns, function(f) {
    tbl <- read.delim(f, row=1); 
    len <- tbl$length;
    cnt <- tbl$expected_count;
    names(len) <- names(cnt) <- rownames(tbl); 
    list(length=len, count=cnt); 
  }); 
  
  ids <- Reduce('union', lapply(d, function(d) names(d[[1]])));
  len <- lapply(d, function(d) d$length); 
  len <- do.call('c', len);
  len <- len[ids]; 
  
  cnm <- TrimPath(fns); 
  cnm <- sub('.tsv', '', cnm); 
  cnt <- matrix(0, nr=length(ids), nc=length(cnm), dimnames=list(ids, cnm)); 
  
  num <- sapply(d, function(d) length(d[[1]])); 
  
  for (i in 1:length(cnm)) cnt[names(d[[i]][[2]]), i] <- as.vector(d[[i]][[2]]); 
  
  list(length=len, number=num, count=cnt); 
}

