### Functions that process BioSystems annotation info and map BioSystems IDs to other IDs

##########################################################################################################################
##########################################################################################################################
# Map gene, protein, etc. to biosystems
Map2Biosystems<-function(bsid, species=c('human'='9606'), to=c('gene', 'protein', 'compound', 'substance', 'pubmed'), by.source=TRUE, 
                         from=c('r', 'gz', 'url'), path=paste(RCHIVE_HOME, 'data/gene.set/public/biosystems', sep='/')) {
  #bsid       Full annotation table of Biosystems downloaded from NCBI website
  #species    One or multiple species to do the mapping; the value is the NCBI Taxonomy ID and the name is the prefix of output files
  #to         Type(s) of entities the Biosystems will be mapped to
  #by.source  Split mapping results by sources of Biosystems if TRUE, such as GO and KEGG
  #from       Source of input data, from previously save R object, downloaded .gz file, or original URL
  #path       Path to the output input files. <path>/r is for all the R objects and <path>/src is for downloaded .gz file
  
  if (!file.exists(path)) dir.create(path);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  
  # all possible input sources
  r<-paste(path, '/r/biosystem2', to, '_fulltable.rds', sep='');
  fnm<-c("biosystems_gene_all.gz", "biosystems_protein.gz", "biosystems_pccompound.gz", "biosystems_pcsubstance.gz", "biosystems_pubmed.gz")
  gz<-paste(path, '/src/', fnm, sep='');
  url<-paste("ftp://ftp.ncbi.nih.gov/pub/biosystems/CURRENT", fnm, sep='/');
  cnm<-c('Gene', 'Protein', 'Compound', 'Substance', 'PubMed');
  names(r)<-names(gz)<-names(url)<-names(cnm)<-c('gene', 'protein', 'compound', 'substance', 'pubmed');
 
  to<-to[to %in% names(r)];
  if (length(to) == 0) to<-'gene';
  r<-r[to];
  gz<-gz[to];
  url<-url[to];
  cnm<-cnm[to];
  
  # Load full mapping data
  frm<-tolower(from)[1];
  ttl<-lapply(to, function(nm) {
    if (frm=='r' & file.exists(r[nm])) readRDS(r[nm]) else {
      if (frm!='gz' | !file.exists(gz[nm])) download.file(url[nm], gz[nm]); # Download data from source
      mp<-read.table(gz[nm], sep='\t', stringsAsFactors=FALSE);
      colnames(mp)<-c('BioSystem_ID', paste(cnm[nm], 'ID', sep='_'), 'Score');
      mp[[1]]<-as.character(mp[[1]]);
      mp[[2]]<-as.character(mp[[2]]);
      saveRDS(mp, file=paste(path, '/r/', 'biosystem2', nm, '_fulltable.rds', sep=''));
      saveRDS(split(mp[[2]], mp[[1]]), file=paste(path, '/r/', 'biosystem2', nm, '_list.rds', sep=''));
      mp;
    }
  });
  names(ttl)<-to;

  sp<-species[species %in% bsid$Taxonomy];
  if (is.null(names(sp))) names(sp)<-sp else names(sp)[names(sp)=='']<-sp[names(sp)==''];
  
  fn<-lapply(names(sp), function(nm) sapply(names(ttl), function(tp) {
    cat('Mapping', tp, 'of', nm, '\n');
    bs<-bsid[bsid$Taxonomy==sp[nm], ];
    t1<-ttl[[tp]];
    t1<-t1[t1[[1]] %in% rownames(bs), ];
    mp<-split(t1[[2]], t1[[1]]);
    
    if (!by.source) saveRDS(mp, file=paste(path, '/r/biosystem2', tp, '_', nm, '.rds', sep='')) else {
      mp0<-split(mp, bsid[names(mp), 1]);
      sapply(names(mp0), function(sc) saveRDS(mp0[[sc]], file=paste(path, '/r/biosystem2', tp, '_', nm, '_',  gsub(' ', '-', sc), '.rds', sep='')))
    }
  }));
  
}

##########################################################################################################################
##########################################################################################################################
# Download and parse general Biosystems annotation information
ParseBiosystemsGeneral<-function(species=c('human'='9606'), ver="ftp://ftp.ncbi.nih.gov/pub/biosystems/CURRENT", download.new=FALSE, 
                          path=paste(RCHIVE_HOME, 'data/gene.set/public/biosystems', sep='/')) {
  # species         Named character vector of NCBI taxanomy ID; the name will be used as prefix of output file
  # ver             The version of BioSystems to download
  # download.new    Whether to re-download source files
  # path            Path to output files
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));
  
  sp<-names(species);
  names(sp)<-species;
  
  ####################################################################################################################
  # Source file names
  fn<-c("biosystems_biosystems_conserved.gz", "biosystems_biosystems_linked.gz", "biosystems_biosystems_similar.gz", 
        "biosystems_biosystems_specific.gz", "biosystems_biosystems_sub.gz", "biosystems_biosystems_super.gz", 
        "biosystems_cdd_specific.gz", "biosystems_gene.gz", "biosystems_gene_all.gz", "biosystems_pccompound.gz", 
        "biosystems_pcsubstance.gz", "biosystems_protein.gz", "biosystems_pubmed.gz", "biosystems_taxonomy.gz", 
        "bsid2info.gz");
  
  # Download source files
  
  if (download.new) download.file('ftp://ftp.ncbi.nih.gov/pub/biosystems/README.txt', paste(path, '/src/README.txt', sep='')); # original README file
  fn0<-paste(ver, fn, sep='/'); # source files
  fn1<-paste(path, 'src', fn, sep='/'); # local files
  dnld<-lapply(1:length(fn), function(i) if (download.new | !file.exists(fn1[i])) download.file(fn0[i], fn1[i]));
 
  #####################
  ### Start parsing ###
  #####################
  
  ###########################################################################################
  # Meta table of bio-systems, row names are the unique BioSystems IDs
  lines<-scan(fn1[grep('bsid2info.gz$', fn1)][1], what='', sep='\n', flush=TRUE);
  split<-strsplit(lines, '\t', useBytes=TRUE);
  bsid<-t(sapply(split, function(s) s[1:8]));
  rownames(bsid)<-bsid[,1];
  bsid[is.na(bsid)]<-'';
  bsid<-data.frame(bsid[, -1], stringsAsFactors=FALSE);
  names(bsid)<-c('Source', 'Accession', 'Name', 'Type', 'Scope', 'Taxonomy', 'Description');
  saveRDS(bsid, file=paste(path, 'r', 'biosystem_all.rds', sep='/'));
  
  # biosystems by source
  cls<-split(bsid[, -1], bsid[, 1]);
  sapply(names(cls), function(nm) saveRDS(cls[[nm]], file=paste(path, '/r/biosystem_', gsub(' ', '-', nm), '.rds', sep='')));
  saveRDS(cls, file=paste(path, 'r', 'biosystem_by_source.rds', sep='/')); 
  
  # Save subset of source-species of selected species (human, mouse, rat, ...)
  tbl<-xtabs(~bsid$Source + bsid$Taxonomy);
  tbl<-tbl[, colnames(tbl) %in% names(sp), drop=FALSE];
  tbl<-tbl[rowSums(tbl)>0, , drop=FALSE];
  sapply(1:nrow(tbl), function(i) sapply(1:ncol(tbl), function(j) {
    fn<-paste(path, '/r/biosystem_', sp[colnames(tbl)[j]], '_', gsub(' ', '-', rownames(tbl)[i]), '.rds', sep='');
    tb<-bsid[bsid$Source == rownames(tbl)[i] & bsid$Taxonomy == colnames(tbl)[j], -c(1, 6), drop=FALSE];
    #file.remove(fn);
    if (nrow(tb) > 0) saveRDS(tb, file=fn);
  }))->nll;
  
  # map organism specific biosystems to corresponding conserved biosystems
  cons2spec<-read.table(fn1[grep("biosystems_biosystems_conserved.gz$", fn1)][1], sep='\t', stringsAsFactors=FALSE);
  saveRDS(split(as.character(cons2spec[,1]), cons2spec[,2]), file=paste(path, 'r', 'biosystem_conserved2specific.rds', sep='/'));
  saveRDS(split(as.character(cons2spec[,2]), cons2spec[,1]), file=paste(path, 'r', 'biosystem_specific2conserved.rds', sep='/'));
  
  cons<-bsid[bsid$Scope=='conserved biosystem', , drop=FALSE];
  cons.id<-split(rownames(cons), cons$Source);
  
  ###########################################################################################
  # Full tables of Biosystems to other ID mapping
  spec2taxo<-read.table(fn1[grep("biosystems_taxonomy.gz$", fn1)][1], sep='\t', stringsAsFactors=FALSE);
  id<-as.character(spec2taxo[[1]]);  
  cons<-as.character(cons2spec[[2]]);
  names(cons)<-cons2spec[[1]];
  spec<-data.frame(Organism=as.character(spec2taxo[[2]]), Conserved=cons[id], Score=spec2taxo[[3]], row.names=id, stringsAsFactors=FALSE);
  spec[is.na(spec)]<-'';
  saveRDS(spec, file=paste(path, 'r', 'biosystem_organism_specific.rds', sep='/'));
  
  tp<-c('gene_all', 'protein', 'pubmed', 'pccompound', 'pcsubstance');
  mp.fn<-sapply(tp, function(tp) fn1[grep(tp, fn1)][1]);
  bs2oth<-lapply(mp.fn, function(fn) read.table(fn, sep='\t', stringsAsFactors=FALSE));
  cnm<-c('Gene', 'Protein', 'PubMed', 'Compound', 'Substance');
  for (i in 1:length(bs2oth)) {
    colnames(bs2oth[[i]])<-c('BioSystem_ID', paste(cnm[i], 'ID', sep='_'), 'Score');
    bs2oth[[i]][[1]]<-as.character(bs2oth[[i]][[1]]);
    bs2oth[[i]][[2]]<-as.character(bs2oth[[i]][[2]]);
    mp<-split(bs2oth[[i]][[2]], bs2oth[[i]][[1]])
    saveRDS(bs2oth[[i]], file=paste(path, '/r/', 'biosystem2', tolower(cnm)[i], '_fulltable.rds', sep=''));
    saveRDS(mp, file=paste(path, '/r/', 'biosystem2', tolower(cnm)[i], '_list.rds', sep=''));
    
    sapply(names(cons.id), function(sc) {
      mp<-mp[names(mp) %in% cons.id[[sc]]];
      saveRDS(mp, file=paste(path, '/r/', 'biosystem2', tolower(cnm)[i], '_conserved_', sc, '.rds', sep=''))
    })
  }
  
  ###########################################################################################
  # Save Log
  # existing full taxonomy table 
  if (download.new) {
    if (file.exists(paste(path, 'r', 'biosystem_all.rds', sep='/'))) {
      bsid0<-readRDS(paste(path, 'r', 'biosystem_all.rds', sep='/'));
      log.old<-list(id=rownames(bsid0), acc=unique(bsid0$Accession), type=unique(bsid0$Type), src=unique(bsid0$Source), taxonomy.old=unique(bsid0$Taxonomy));        
    } else {
      log.old<-list(id=c(), acc=c(), type=c(), src=c(), taxonomy=c())
    }
    log.new<-list(id=rownames(bsid), acc=unique(bsid$Accession), type=unique(bsid$Type), src=unique(bsid$Source), taxonomy.old=unique(bsid$Taxonomy));
    names(log.old)<-names(log.new)<-c('ID', 'Accession', 'Type', 'Source', 'Taxonomy');
    
    # updates
    up<-list(
      N = sapply(log.new, length),
      Added = lapply(1:5, function(i) setdiff(log.new[[i]], log.old[[i]])),
      Removed = lapply(1:5, function(i) setdiff(log.old[[i]], log.new[[i]]))
    )
    
    # update logs
    log<-readRDS(paste(path, 'log.rds', sep='/'));
    log<-c(log, list(up));
    names(log)[length(log)]<-as.character(Sys.Date());
    saveRDS(log, file=paste(path, 'log.rds', sep='/'));
  }
}

