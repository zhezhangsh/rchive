# Parse the NCBI Entrez genes of multiple species
ParseEntrez<-function(species=c('human'='9606'), download.new=TRUE, path=paste(RCHIVE_HOME, 'data/gene/public/entrez', sep='/')) {
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));
 
  fn<-paste(path, 'r', 'all_genes.rds', sep='/');
  
  # Re-download and load all gene table
  if (download.new | !file.exists(fn)) {
    gz<-paste(path, 'src', 'All_Data.gene_info.gz', sep='/');
    download.file("ftp://ftp.ncbi.nih.gov/gene//DATA/GENE_INFO/All_Data.gene_info.gz", gz);
    
    hd<-scan(gz, nline=1, what='', sep='\n', flush=TRUE);
    hd<-strsplit(hd, ' ')[[1]][2:16];
    
    ### Give error: could not allocate memory (2048 Mb) in C function 'R_AllocStringBuffer'
    #all<-read.table(fn, sep='\t', head=FALSE, quote='\"', comment='', skip=1, stringsAsFactors=FALSE);
    
    gn<-scan(fn, sep='\n', what='', flush=TRUE, skip=1);
    gn<-strsplit(gn, '\t');
    all<-do.call('rbind', gn);
    all<-data.frame(all, stringsAsFactors=FALSE);
    names(all)<-hd;
    
    saveRDS(all, file=fn);
  } else all<-readRDS(fn);
  
  ############################################################
  # Parse the NCBI Entrez genes of a species
  parser<-function(all, tax, prefix, path) {
    # all     The full Entrez gene table
    # tax     Taxonomy ID
    # prefix  Species name as file prefix
    # path    Path to output files
        
    gn<-all[all[[1]]==tax, , drop=FALSE];
    
    if (nrow(gn) > 0) {
      rownames(gn)<-as.vector(gn$GeneID);
      gn<-gn[, !(names(gn) %in% c('GeneID', 'tax_id'))];
      
      # Map Entrez ID to official symbol
      sym<-as.vector(gn$Symbol);
      names(sym)<-rownames(gn);
      sym<-sym[!is.na(sym) & sym!='-'];
      id2sym<-sym;
      sym2id<-names(sym);
      names(sym2id)<-sym;
      
      # Map Entrez ID to symbol synonyms
      syn<-as.matrix(gn[, colnames(gn) %in% c('Symbol', 'Synonyms', 'Symbol_from_nomenclature_authority')]);
      syn<-lapply(1:nrow(syn), function(i) unlist(strsplit(syn[i,], '\\|'), use.names=FALSE));
      names(syn)<-rownames(gn);
      id2syn<-lapply(syn, function(x) unique(x[x!='-']));
      syn2id<-split(rep(names(id2syn), sapply(id2syn, length)), unlist(id2syn, use.names=FALSE));
      syn2id<-lapply(syn2id, unique);
      
      # Formatted DataTable for Shiny rendering
      gn.tbl<-data.frame(
        ID = paste('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=', rownames(gn), '" target="_blank">', rownames(gn), '</a>', sep=''),
        Symbol = as.vector(gn$Symbol),
        Location = as.vector(gn$map_location),
        Type = as.vector(gn$type_of_gene),
        Description = as.vector(gn$description)
      )
      
      # save outputs
      saveRDS(gn, file=paste(path, '/r/', prefix, '_genes_full.rds', sep=''));
      saveRDS(gn.tbl, file=paste(path, '/r/', prefix, '_genes_datatable.rds', sep=''));
      saveRDS(id2sym, file=paste(path, '/r/', prefix, '_genes_id2symbol.rds', sep=''));
      saveRDS(sym2id, file=paste(path, '/r/', prefix, '_genes_symbol2id.rds', sep=''));
      saveRDS(id2syn, file=paste(path, '/r/', prefix, '_genes_id2synonyms.rds', sep=''));
      saveRDS(syn2id, file=paste(path, '/r/', prefix, '_genes_synonyms2id.rds', sep=''));
    }
    
    rownames(gn);
  }
  ############################################################

  ids<-lapply(names(species), function(nm) parser(all, species[nm], nm, path));
  if (download.new) {
    log<-readRDS(paste(path, 'log.rds', sep='/'));
    log<-c(log, list(up));
    names(log)[length(log)]<-as.character(Sys.Date());
    saveRDS(log, paste(path, 'log.rds', sep='/'));
  }
}
