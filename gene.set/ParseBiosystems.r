ParseBiosystems<-function(species=c('human'='9606'), ver="ftp://ftp.ncbi.nih.gov/pub/biosystems/CURRENT", download.new=FALSE, 
                          path=paste(RCHIVE_HOME, 'data/gene.set/public/biosystems', sep='/')) {
  # species         Named character vector of NCBI taxanomy ID; the name will be used as prefix of output file
  # ver             The version of BioSystems to download
  # download.new    Whether to re-download source files###
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
  download.file('ftp://ftp.ncbi.nih.gov/pub/biosystems/README.txt', paste(path, '/src/README.txt', sep='')); # original README file
  fn0<-paste(ver, fn, sep='/'); # source files
  fn1<-paste(path, 'src', fn, sep='/'); # local files
  dnld<-lapply(1:length(fn), function(i) if (download.new | !file.exists(fn1[i])) download.file(fn0[i], fn1[i]));
 
  ####################################################################################################################
  # Meta table of bio-systems, row names are the unique BioSystems IDs
  lines<-scan(fn1[grep('bsid2info.gz$', fn1)], what='', sep='\n', flush=TRUE);
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
  
  ##################################################################################
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

#############################################################
#############################################################

# Selected species
sp<-c(
  '9606' = 'human',
  '10090' = 'mouse',
  '10116' = 'rat',
  '9598' = 'chimp',
  '9823' = 'pig',
  '9031' = 'chicken',
  '9615' = 'dog',
  '9913' = 'cow',
  '6239' = 'worm',
  '7227' = 'fly',
  '7955' = 'zebrafish',
  '4932' = 'yeast',
  '511145' = 'ecoli'
);





###########################################################################################
# BioSystems to gene mapping
bs2gn<-read.table(paste(path, 'src', 'biosystems_gene_all.gz', sep='/'), sep='\t');
names(bs2gn)<-c('BioSystem_ID', 'Gene_ID', 'Score');
bs2gn[[1]]<-as.character(bs2gn[[1]]);
bs2gn[[2]]<-as.character(bs2gn[[2]]);
saveRDS(bs2gn, file=paste(path, 'r', 'biosystem-gene_all.rds', sep='/'));

# BioSystems to gene as list
b2g<-split(bs2gn[,2], bs2gn[,1]);
saveRDS(b2g, file=paste(path, 'r', 'biosystem2gene_all.rds', sep='/'));

# Split list by source
id<-intersect(names(b2g), rownames(bsid));
b2g<-b2g[id];
bsid1<-bsid[id, ];
cls.gn<-split(b2g, bsid1$Source);
sapply(names(cls.gn), function(nm) saveRDS(cls.gn[[nm]], file=paste(path, '/r/biosystem2gene_', gsub(' ', '-', nm), '.rds', sep='')));

# Save sub-list of source-species of selected species (human, mouse, rat, ...)
sapply(1:nrow(tbl), function(i) sapply(1:ncol(tbl), function(j) {
  fn<-paste(path, '/r/biosystem2gene_', sp[colnames(tbl)[j]], '_', gsub(' ', '-', rownames(tbl)[i]), '.rds', sep='');
  tb<-bsid1[bsid1$Source == rownames(tbl)[i] & bsid1$Taxonomy == colnames(tbl)[j], -c(1, 6), drop=FALSE];
  #file.remove(fn);
  if (nrow(tb) > 0) saveRDS(b2g[rownames(tb)], file=fn);
}))->nll;
