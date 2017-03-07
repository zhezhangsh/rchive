FormatGenesetCollection <- function(name, path, meta, mapp, species, 
                                    path.gene=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez/r', sep='/')) {
  cnm <- c("Collection", "Name", "Species", "Size", "URL");
  if (!identical(cnm, colnames(meta))) stop('Error: invalid metadata column(s)\n'); 
  if (!identical(names(mapp), rownames(meta))) stop('Error: IDs not match, metadata vs. lists\n');
  if (length(mapp) != length(names(mapp))) stop('Error: Duplicated IDs\n');

  require(ff); 
  
  path <- paste(path, name, sep='/'); 
  if (!dir.exists(path)) dir.create(path, recursive = TRUE); 
  
  # files to gene identifier mapping
  fs <- paste(path.gene, paste(species, '_genes_id2symbol.rds', sep=''), sep='/');
  fe <- paste(path.gene, paste(species, '_genes_id2ensg.rds', sep=''), sep='/');
  names(fs) <- names(fe) <- species; 
  
  # Format and save metadata
  f1 <- paste(path, 'metadata.txt', sep='/');
  f2 <- paste(path, 'metadata.ff', sep='/');
  f3 <- paste(path, 'metadata.rds', sep='/');
  write.csv2(cbind(ID=rownames(meta), meta), file=f1, quote = TRUE, row.names = FALSE); 
  ff <- read.csv2(f1, row=1, header = TRUE, stringsAsFactors = TRUE);
  saveRDS(ff, f2); 
  saveRDS(meta, f3); 
  
  # Split lists by colleciton and species
  cll <- split(mapp, meta$Collection); 
  smm <- lapply(names(cll), function(nm1) {
    mp <- cll[[nm1]]; 
    mt <- meta[names(mp), , drop=FALSE]; 
    sp <- split(mp, mt$Species); 
    sp <- sp[names(sp) %in% species]; 
    pt <- paste(path, nm1, sep='/'); 
    if (!dir.exists(pt)) dir.create(pt, recursive = TRUE); 
    sm <- lapply(names(sp), function(nm2) {
      m <- sp[[nm2]]; 
      p <- paste(path, nm1, nm2, sep='/'); 
      if (!dir.exists(p)) dir.create(p, recursive = TRUE); 
      cat(nm1, nm2, length(m), '\n'); 
      d <- names(m); 
      x <- rep(names(m), sapply(m, length)); 
      y <- unlist(m, use.names=FALSE); 
      i <- which(x!=''&y!=''&!is.na(x)&!is.na(y));
      x <- x[i]; 
      y <- y[i]; 
      m <- split(y, x); 
      m <- lapply(m, unique); 
      m <- lapply(m, function(m) m[order(as.numeric(m))]); 
      m <- m[d];
      s <- 
      names(m) <- d; 
      saveRDS(m, paste(p, 'entrez.rds', sep='/')); 
      saveRDS(meta[meta$Species==nm2&meta$Collection==nm1, , drop=FALSE], paste(p, 'metadata.rds', sep='/')); 
      
      n1 <- sapply(m, length);
      n2 <- rep(0, length(n1));
      n3 <- rep(0, length(n1));
      bg <- unique(as.vector(y));
      bg <- bg[order(as.numeric(bg))]; 
      bg <- list(entrez=bg); 
      
      # map to official gene symbole
      if (file.exists(fs[nm2])) {
        cat('Mapping to gene symbol', '\n');
        s <- readRDS(fs[nm2]); 
        z <- s[y]; 
        u <- unique(as.vector(z));
        z <- split(z[!is.na(z)], x[!is.na(z)]);
        z <- lapply(z, unique); 
        z <- lapply(z, sort); 
        z <- z[d]; 
        names(z) <- d; 
        n2 <- sapply(z, length); 
        bg$symbol <- sort(u[u!=''&!is.na(u)]); 
        saveRDS(z, paste(p, 'symbol.rds', sep='/')); 
      };
      # map to ensembl gene ID
      if (file.exists(fe[nm2])) {
        cat('Mapping to ensembl id', '\n');
        s <- readRDS(fe[nm2]); 
        z <- s[y]; 
        z <- split(z, x);
        z <- lapply(z, function(z) unlist(z, use.names=FALSE)); 
        z <- lapply(z, function(z) if (!is.null(z)) z[!is.na(z)] else NULL); 
        z <- lapply(z, unique); 
        z <- lapply(z, function(z) if(!is.null(z)) sort(z) else NULL); 
        z <- z[d]; 
        names(z) <- d; 
        n3 <- sapply(z, length); 
        bg$ensembl <- sort(unique(unlist(z, use.names=FALSE)));  
        saveRDS(z, paste(p, 'ensembl.rds', sep='/')); 
      };
      ns <- cbind(Num_Entrez=n1, Num_Symbol=n2, Num_Ensembl=n3); 
      rownames(ns) <- d;
      saveRDS(bg, paste(p, 'background.rds', sep='/')); 
      saveRDS(ns, paste(p, 'count.rds', sep='/')); 
      list(background=bg, count=ns); 
    });
    names(sm) <- names(sp); 
    saveRDS(sm, file=paste(path, nm1, 'summary.rds', sep='/')); 
    sm; 
  }); 
  names(smm) <- names(cll); 
  
  saveRDS(smm, file=paste(path, 'summary.rds', sep='/')); 
  
  smm; 
}