# Get details of Entrez genes from NCBI website by parsing from xml files

####################### THIS IS A SLOW PROCESS, REQUIRES DOWNLOADING SEVERAL GB OF DATA FROM NCBI ################################
#  Dependent on: saved gene annotation downloaded from NCBI FTP server; an R object named "anno" whose row names are Entrez IDs

GetEntrezDetail<-function(species='human', num.clusters=1, block.size=500, path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene/public/entrez', sep='/')) {
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
  if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);
  
  ########################################################################################################################
  ########################################################################################################################
  ## Parse a set of NCBI Entrez genes
  parse.set<-function(id) {
    ########################################################################################################################
    ## clean up HTML tags in a line
    cleanHtml<-function(l) {
      l<-gsub("<.*?>", " ", l); # remove open tags
      l<-gsub("</.*?>", " ", l); # remove close tags
      l<-gsub(" +", " ", l); # replace multiple spaces with one
      l<-sub('^ ', '', l); # remove space at the begining of the line
      l<-sub(' $', '', l); # remove space at the end of the line
      l;
    }
    ########################################################################################################################librar
    
    ########################################################################################################################
    ## Parse one NCBI Entrez gene
    parse.one<-function(xml) {
      ####################################################################################
      # getting gene location of both genome versions
      ind.ver<-grep('<Gene-commentary_heading>GRCh3[78]</Gene-commentary_heading>', xml); # lines with assembly version (hg19 or hg20)
      ind.open<-grep('<Gene-commentary>', xml);
      ind.close<-grep('</Gene-commentary>', xml);
      
      # Get blocks of xml lines with location information
      blocks<-lapply(ind.ver, function(ind) {
        ind1<-ind.open[ind.open>ind];
        ind2<-ind.close[ind.close>ind];
        if (length(ind1)>0 & length(ind2)>0) {
          if (min(ind1) < min(ind2)) xml[ind:min(ind2)] else NA
        } else NA;
      });
      blocks<-blocks[sapply(blocks, length)>1];
      
      # Get location in each block
      locs<-lapply(blocks, function(b) {
        # get gene strand
        str<-b[grep('<Seq-interval_strand>', b)+1];
        str<-sapply(str, function(str) if (grepl("plus", str)) str<-'+' else if (grepl("minus", str)) str<-'-' else str<-'*');
        
        # parse location info
        ln<-list(
          Version=b[grep('<Gene-commentary_heading>GRCh3', b)],
          Chromosome=b[grep('<Gene-commentary_label>Chr', b)],
          From=b[grep('<Seq-interval_from>', b)],
          To=b[grep('<Seq-interval_to>', b)]
        );
        
        lc<-lapply(ln, cleanHtml); # clean HTML tags and extra spaces
        
        # re-formatting
        lc$Chromosome<-sub('Chromosome ', '', lc$Chromosome);
        lc$From<-as.integer(lc$From);
        lc$To<-as.integer(lc$To);
        lc$Strand<-as.vector(str);
        
        lc;
      });
      names(locs)<-sapply(locs, function(x) x$Version);
      ####################################################################################
      
      list(
        "ID"=cleanHtml(xml[grep('<Gene-track_geneid>', xml)]), # gene ID
        "Official Symbol"=cleanHtml(xml[grep('<Gene-ref_locus>', xml)]), # official symbol
        "Synonyms"=cleanHtml(xml[grep('<Gene-ref_syn_E>', xml)]), # Synonyms
        "Name"=cleanHtml(xml[grep('<Gene-ref_desc>', xml)]), # full name
        "Cytoband"=cleanHtml(xml[grep('<Gene-ref_maploc>', xml)]), # cytoband
        "Location"=locs,
        "Summary"=cleanHtml(xml[grep('<Entrezgene_summary>', xml)]) # summary
      );
    }
    ########################################################################################################################
    
    library(RCurl)
    
    url<-paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&rettype=xml&id=", paste(id, collapse=','), sep='');
    xmls<-getURL(url); # downlaod data from NCBI
    xmls<-strsplit(xmls, '\n')[[1]]; # break lines
    
    # split xml file by genes
    ind1<-grep('<Entrezgene>', xmls);
    ind2<-grep('</Entrezgene>', xmls);
    if (length(ind1)!=length(ind2) | length(which(ind2<=ind1))>0 | length(which(ind1[-1]<=ind2[-length(ind2)]))>0)
      stop('Error: spliting xml by genes.\n');
    xmls<-lapply(1:length(ind1), function(i) xmls[ind1[i]:ind2[i]]);
    #names(xmls)<-id;
    
    ######################################
    lapply(xmls, parse.one); # parsing
    ######################################
  }
  ########################################################################################################################
  ########################################################################################################################
  
  # Get all Entrez IDs from previous processed gene annotation based downloaded data from NCBI FTP server
  fn.anno<-paste(path, '/r/', species, '_genes_full.rds', sep='');   # Required file
  if (!file.exists(fn.anno)) stop("Error: annotation file for specie", species, "not exists.\n");
  anno<-readRDS(fn.anno);
  all.ids<-rownames(anno); # Entrez gene ID
  
  fn.old <- paste(path, '/r/', species, '_by_id.rds', sep='');
  if (file.exists(fn.old)) gn.old <- readRDS(fn.old) else gn.old <- c();
  id.new <- rownames(anno)[!(rownames(anno) %in% names(gn.old))];
  
  ids <- split(id.new, rep(1:ceiling(length(id.new)/block.size), each=block.size)[1:length(id.new)]); # split into smaller blocks
  
  # fetch and parse gene info in batch
  if (num.clusters > 1) { # using parallel processing
    library(snow);
    cl<-makeCluster(4, type='SOCK');
    gn<-clusterApplyLB(cl, ids, parse.set);
    stopCluster(cl);
  } else {
    gn<-list();
    for (i in 1:length(ids)) {
      print(i);
      gn[[i]]<-parse.set(ids[[i]]);
    }
  }
  
  # merge of 
  gn.all<-do.call('c', gn);
  names(gn.all)<-sapply(gn.all, function(g) g$ID);
  
  # insert URL
  gn.all<-lapply(gn.all, function(gn) {
    id<-gn[[1]];
    url<-c('URL'=paste('http://www.ncbi.nlm.nih.gov/gene/?term=', id, sep=''));
    append(gn, url, after=4);
  });
  gn.all <- c(gn.old, gn.all); 
  gn.all<-gn.all[all.ids];
  
  saveRDS(gn.all, file=paste(path, '/r/', species, '_by_id.rds', sep=''));
}