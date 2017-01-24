# Map probe sequences to reference transcriptome
MapProbe2Gene <- function(seq, bsgenome, txdb) {
  # seq       A set of the DNA sequences of probes, all having the same length
  # bsgenome  BSgenome object of reference genome
  # txdb      TxDb object of reference transcriptome
  
  require(Biostrings); 
  require(BSgenome); 
  require(GenomicFeatures);
  
  names(seq) <- 1:length(seq); 
  nc  <- nchar(seq); 
  seq <- seq[nc==max(nc)]; 
  seq <- DNAStringSet(seq); 
  pd  <- PDict(seq); 
  
  tx <- extractTranscriptSeqs(bsgenome, txdb);
  
  mp <- lapply(tx, function(tx) {
    m <- matchPDict(pd, tx); 
    names(m)[elementLengths(m)>0]; 
  }); 
  mp <- mp[elementLengths(mp)>0]; 
  
  txs <- transcripts(txdb, columns=c('tx_name', 'gene_id')); 
  txs <- txs[txs$tx_name %in% names(mp)]; 
  tnm <- as.vector(txs$tx_name); 
  gid <- as.list(txs$gene_id); 
  tnm <- rep(tnm, elementLengths(gid)); 
  mps <- mp[tnm]; 
  names(mps) <- unlist(gid); 
  mps <- split(mps, names(mps)); 
  mps <- lapply(mps, function(x) unique(unlist(x, use.names=FALSE))); 
  
  seq2gn <- rep(names(mps), elementLengths(mps)); 
  names(seq2gn) <- seq[as.integer(unlist(mps, use.names=FALSE))]; 
  seq2gn <- split(as.vector(seq2gn), names(seq2gn)); 
  
  seq2gn; 
}