##############################################
# Mapping between KEGG databases and entries #
##############################################

MapKegg2Gene.2 <- function(types=, species=c('human'='hsa'), 
                           path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg', sep='/')) {
  require(RCurl); 
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  
  nm<-names(species);
  if (is.null(nm)) names(species)<-species else names(species)[is.na(names(species))]<-species[is.na(names(species))];
  
  # ID types
  types<-tolower(types);

  grp <- c('compound'="cpd", 'dgroup'="dg", 'drug'="dr", 'disease'="ds", 'enzyme'="ec",
           'module'="md", 'omim'="omim", 'pathway'="path", 'pfam'="pf", 'reaction'="rn"); 
  
  ################################################################
  # annotation
  url <- paste('http://rest.kegg.jp/list', names(grp), sep='/');
  names(url) <- names(grp);
  lapply(names(url), function(nm) {
    if (url.exists(url[nm])) {
      l <- readLines(url[nm]); 
      l <- strsplit(l, '\t');
      x <- sapply(l, function(x) x[1]);
      y <- sapply(l, function(x) x[2]);
      x <- sapply(strsplit(x, ':'), function(x) x[2]); 
      names(y) <- x;
      df <- data.frame(Name=y, URL=paste('http://www.kegg.jp/dbget-bin/www_bget?', x, sep=''), stringsAsFactors = FALSE);
      rownames(df) <- x;
      fnm <- paste(path, '/r/anno_', nm, '.rds', sep='');
      saveRDS(df, fnm); 
      fnm;
    } else '';
  }) -> x;
  ################################################################
  ################################################################
  # Genes
  url <- paste('http://rest.kegg.jp/list', species, sep='/');
  ################################################################
  
  if (length(types)>0) grp <- grp[names(grp) %in% types];
  
  ids <- lapply(names(species), function(nm) {
    cat(nm, '\n'); 
    spe <- species[nm]; 
    gns <- names(ConvertGene2EntrezKeggApi(spe));
    gns <- split(gns, rep(1:ceiling(length(gns)/100), each=100)[1:length(gns)]);
    
    mps <- lapply(gns, function(g) {
      l <- readLines(paste('http://rest.genome.jp/link', paste(g, collapse='+'), sep='/')); 
      l <- strsplit(l, '\t'); 
      x <- sapply(l, function(x) x[2]);
      x <- strsplit(x, ':'); 
      y <- sapply(x, function(x) x[2]); 
      z <- sapply(x, function(x) x[1]); 
      names(y) <- sub(paste(spe, ':', sep=''), '', sapply(l, function(l) l[1]));
      map <- split(y, z); 
      map <- map[names(map) %in% grp];
      names(map) <- sapply(names(map), function(x) names(grp)[grp==x]);
      map; 
    }); 
    
    cll <- lapply(mps, names);
    cll <- sort(unique(unlist(cll)));
    ids <- lapply(cll, function(c) {
      cat(c, '\n'); 
      ttl <- lapply(mps, function(m) if (c %in% names(m)) m[[c]] else NULL);
      gns <- unlist(lapply(ttl, names), use.names=FALSE);
      ids <- unlist(lapply(ttl, as.vector), use.names=FALSE);
      map <- split(gns, ids);
      map <- lapply(map, unique);
      fnm <- paste(path, '/r/', nm, '_', c, '2gene.rds', sep='');
      saveRDS(map, fnm); 
      cat(fnm, '\n'); 
      names(map); 
    }); 
    
    names(ids) <- cll;
    ids;
  });
  
  names(ids) <- names(species);
  
  lapply(names(ids), function(spe) {
    id <- ids[[spe]]; 
    nm <- lapply(names(id), function(typ) {
      i <- id[[typ]]; 
      u <- paste()
    }); 
  });
  
  ids;
}

######################################################################################################
######################################################################################################
# Map types of KEGG IDs to Entrez genes of one or more species
MapKegg2Gene<-function(types=c('pathway'), species=c('human'='hsa'), path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg', sep='/')) {
  # species   Named character vector of species codes; the name will be used as prefix of output file
  # path      Path to output files
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  
  nm<-names(species);
  if (is.null(nm)) names(species)<-species else names(species)[is.na(names(species))]<-species[is.na(names(species))];
  
  # ID types
  dbs<-strsplit("pathway | brite | module | ko | compound | glycan | reaction | rpair | rclass | enzyme | disease | drug | dgroup | environ", ' \\| ')[[1]];
  types<-tolower(types);
  types<-types[types %in% dbs];
  
  map2entrez<-lapply(species, ConvertGene2EntrezKeggApi);
  
  if (length(types)==0 | length(species)==0) c() else {
    ids<-lapply(types, function(tp) {
      ids<-lapply(names(species), function(nm) {
        cat(tp, '\t', nm, '\n');
        sp<-species[nm];
        
        # Full names
        if (tp=='pathway' | tp=='module' ) sp0<-sp else sp0<-'';
        anno<-ListEntryKeggApi(tp, sp0);
        id<-sapply(strsplit(rownames(anno), ':'), function(x) x[length(x)]);
        desc<-sapply(strsplit(anno[[1]], ' - '), function(x) x[1]);
        names(desc)<-id;
        
        # Mapping
        mapped<-LinkEntryKeggApi(sp, tp);
        if (nrow(mapped) > 0) {
          id0<-sapply(strsplit(mapped[[1]], ':'), function(x) x[length(x)]);
          mp<-split(mapped[[2]], id0);
          mp<-lapply(mp, function(x) {
            gn<-unlist(map2entrez[[nm]][x], use.names=FALSE);
            if (length(gn) == 0) c() else unique(gn[!is.na(gn)])
          });
          mp<-mp[sapply(mp, length)>0];
          desc<-desc[names(desc) %in% names(mp)];
          mp<-mp[names(desc)];
          
          # Save results
          fn<-paste(path, '/r/', nm, '_', tp, '2gene.rds', sep='');
          saveRDS(list(Organism=species[nm], Name=desc, Map=mp), file=fn);
          
          id;
        } else c();  
      });
      names(ids)<-names(species);
      ids;
    });
    names(ids)<-types;
    ids;
  }
}

######################################################################################################
######################################################################################################
# Usee the KEGG "ID conversion" API to convert all KEGG gene ID of a species to NCBI gene ID
ConvertGene2EntrezKeggApi<-function(org) {
  # org: Organism code, such as 'hsa' or 'cel', for the alternative query form as above
  url<-paste("http://rest.kegg.jp/conv", org, "ncbi-geneid", sep='/');
  cat(url, '\n');
  
  download.file(url, 'tmp.csv');
  mp<-read.csv2('tmp.csv', sep='\t', header=FALSE, stringsAsFactors=FALSE);
  id<-sub("^ncbi-geneid:", '', mp[,1]);
  split(id, mp[,2]);
}

######################################################################################################
######################################################################################################
# Use the KEGG "Entry list" API to get annotation of a database
# URL form: http://rest.kegg.jp/list/<database>
#     <database> = pathway | brite | module | ko | genome | <org> | compound | glycan | reaction | rpair | rclass | enzyme | disease | drug | dgroup | environ | organism
#     <org> = KEGG organism code or T number
# Or: http://rest.kegg.jp/list/<database>/<org>
#     <database> = pathway | module
#     <org> = KEGG organism code
ListEntryKeggApi<-function(db, org='') {
  # db    Database name; see "http://www.kegg.jp/kegg/rest/keggapi.html" for full list
  # org   Organism code, such as 'hsa', for the alternative query form as above
  db<-tolower(db[1]);
  org<-tolower(org[1]);
  
  # Build URL
  url<-paste("http://rest.kegg.jp/list", db, org, sep='/');
  cat(url, '\n');
  
  download.file(url, 'tmp.csv');
  anno<-read.csv2('tmp.csv', sep='\t', header=FALSE, stringsAsFactors=FALSE);
  rownames(anno)<-anno[[1]];
  anno<-anno[, -1, drop=FALSE];
  file.remove('tmp.csv');
  
  anno;
}

######################################################################################################
######################################################################################################
# Use the KEGG "Linked Entry" API to map IDs of 2 databases
# URL form, full databases: http://rest.kegg.jp/link/<target_db>/<source_db>[/<option>]
#     <target_db> = <database>
#     <source_db> = <database>
#     <database> = pathway | brite | module | ko | genome | <org> | compound | glycan | reaction | rpair | rclass | enzyme | disease | drug | dgroup | environ
#     <option> = turtle | n-triple
# URL form, specify entries: http://rest.kegg.jp/link/<target_db>/<dbentries>[/<option>]
#     <dbentries> = KEGG database entries involving the following <database>
#     <database> = pathway | brite | module | ko | genome | <org> | compound | glycan | reaction | rpair | rclass | enzyme | disease | drug | dgroup | environ | genes
#     <option> = turtle | n-triple
LinkEntryKeggApi<-function(to, from, option='') {
  # to      Target database, target_db
  # from    Source database name, or a vector of source database entries
  # option  turtle | n-triple
  
  url<-paste("http://rest.kegg.jp/link", to, paste(from, collapse='+'), option, sep='/');
  
  download.file(url, 'tmp.csv');
  
  if (file.info('tmp.csv')[,1]<=1) mp<-data.frame(matrix(nr=0, nc=2)) else mp<-read.csv2('tmp.csv', sep='\t', header=FALSE, stringsAsFactors=FALSE);
  colnames(mp)<-c('From', 'To');
  file.remove('tmp.csv');
  
  mp;
}

#######################################
#######################################
######### Obsolete functions ##########
#######################################
#######################################


######################################################################################################
######################################################################################################
# Use the KEGG dbget engine to map a ID to IDs of all KEGG database
# Download and parse from an HTML page to get the mapping results
# The URL of the HTML page will be built as:
#     http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+<Abbr>:<ID> 
#         - <Abbr> is the abbrevation of the database defining the query ID
#         - <ID> is the query ID
# For example, "http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+path:hsa01100" maps pathway hsa01100 to IDs of all databases
# Refer to this page for available databases and name abbreviations: http://www.genome.jp/dbget/
MapKeggOne2All<-function(id, type, verbose=FALSE) {
  # id    A unique KEGG ID; gene, compound, drug, etc.
  # type  Type of the ID; such as 'gn', 'cpd', etc. see http://www.genome.jp/dbget
  
  #library(RCurl);
  library(devtools);
  install_github("zhezhangsh/rchive");
  library(rchive);
  
  # URL prefix
  prefix<-"http://www.genome.jp/dbget-bin/get_linkdb?-t+alldb+";
  
  types<-c(
    'pathway' = 'path', 'path' = 'path',
    'brite' = 'br', 'br' = 'br',
    'module' = 'md', 'md' = 'md',
    'enzyme' = 'ec', 'ec' = 'ec',
    'reaction' = 'rn', 'rn' = 'rn',
    'rpair' = 'rp', 'rp' = 'rp', 
    'rclass' = 'rc', 'rc' = 'rc', 
    'orthology' = 'ko', 'ko' = 'ko', 
    'gene' = 'gn', 'gn' = 'gn', 
    'disease' = 'ds', 'ds' = 'ds',
    'compound' = 'cpd', 'cpd' = 'cpd',
    'glycan' = 'gl', 'gl' = 'gl',
    'drug' = 'dr', 'dr' = 'dr',
    'dgroup' ='dg', 'dg' = 'dg',
    'environ' = 'ev', 'ev' = 'ev'
  )
  tp<-types[tolower(type)];
  if (is.na(tp)) tp<-type;
  
  if (verbose) cat(tp, '=', id, '\n');
  
  # Read in HTML file from URL
  url<-paste(prefix, tp, ':', id, sep='');
  #html<-strsplit(getURL(url), '\n')[[1]];
  download.file(url, 'tmp.html');
  html<-scan('tmp.html', sep='\n', what='', flush=TRUE);
  file.remove('tmp.html');
  
  if (verbose) {
    cat(url, '\n');
    cat(length(html), 'lines\n');
  }
  ind.db<-grep('^<span ', html); # lines specifying source databases 
  
  if (length(ind.db) == 0) list() else {
    html<-gsub('</a> ', '\t</a>', html);
    ln<-CleanHtmlTags(html); # Remove HTML tags
    ln<-gsub('\t +', '\t', ln);
    ln<-gsub('\t$', '\t ', ln);
    
    # Mapping
    nm.db<-ln[ind.db];
    n<-c(ind.db[-1]-ind.db[-length(ind.db)], length(ln)-ind.db[length(ind.db)]+1);
    nm<-rep(c('', nm.db), c(ind.db[1]-1, n));
    mp<-split(ln, nm)[-1];
    mp<-lapply(mp, function(mp) mp[-1]);
    mp<-lapply(mp, function(mp) mp[mp!='' & !is.na(mp)]);
    mp<-mp[sapply(mp, length)>0];
    
    # Turn to named vectors
    mp<-lapply(mp, function(mp) {
      x<-do.call('rbind', strsplit(mp, '\t'));
      y<-x[,2];
      names(y)<-x[,1];
      y[y==' ' | is.na(y)]<-'';
      y;
    })
  
    mp;
  }
}


######################################################################################################
######################################################################################################
# Use the KEGG dbget engine to map a ID to IDs of specified KEGG database
# Download and parse from an HTML page to get the mapping results
# The URL of the HTML page will be built as:
#     http://www.genome.jp/dbget-bin/get_linkdb?-t+<DB>+<Abbr>:<ID>
#         - <DB> is the name of the database to be mapped to
#         - <Abbr> is the abbrevation of the database defining the query ID
#         - <ID> is the query ID
# For example, "http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:hsa01100" maps pathway hsa01100 to genes
# Refer to this page for available databases, database names, name abbreviations: http://www.genome.jp/dbget/
MapKeggOne2One<-function(id, type, to, verbose=FALSE) {
  # id    A unique KEGG ID, gene, compound, drug, etc.
  # type  Type of the ID, use database abbreviation such as 'gn', 'path', and 'cpd', see http://www.genome.jp/dbget
  # to    The database that the id will be mapped to, use database name, such as 'genes', 'pathway', and 'compound'
  
  library(devtools);
  install_github("zhezhangsh/rchive");
  library(rchive);
  
  # URL prefix
  url<-paste("http://www.genome.jp/dbget-bin/get_linkdb?-t+", to, '+', type, ':', id, sep='');
  
  download.file(url, 'tmp.html');
  html<-scan('tmp.html', sep='\n', what='', flush=TRUE);
  file.remove('tmp.html');
  
  if (verbose) {
    cat(url, '\n');
    cat(length(html), 'lines\n');
  }
  
  # Parsing HTML file
  ln<-html[(grep('^---', html)+1):(grep('integrated database retrieval system', html)-1)];
  ln<-gsub('</a> ', '\t</a>', ln);
  ln<-CleanHtmlTags(ln);
  ln<-gsub('\t +', '\t', ln);
  ln<-gsub('\t$', '\t ', ln);
  
  # Turn to named vectors
  x<-do.call('rbind', strsplit(ln, '\t'));
  mp<-x[,2];
  names(mp)<-x[,1];
  mp[mp==' ' | is.na(mp)]<-'';
  
  mp;
}


# Map all KEGG pathways to Entrez genes of one or more species
MapKeggPath2Gene<-function(species=c('human'='hsa'), path=paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg', sep='/')) {
  # species   Named character vector of species codes; the name will be used as prefix of output file
  # path      Path to output files
  
  if (!file.exists(path)) dir.create(path, recursive=TRUE);
  if (!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
  
  nm<-names(species);
  if (is.null(nm)) names(species)<-species else names(species)[is.na(names(species))]<-species[is.na(names(species))];
  
  ids<-lapply(names(species), function(nm) {
    cat(nm, '\n');
    sp<-species[nm];
    
    # Pathway full names
    anno<-ListEntryKeggApi('pathway', sp);
    id<-sub('^path:', '', rownames(anno));
    desc<-sapply(strsplit(anno[[1]], ' - '), function(x) x[1]);
    names(desc)<-id;
    
    # Mapping
    mapped<-LinkEntryKeggApi(sp, 'pathway');
    pth<-sub('^path:', '', mapped[[1]]);
    #gn<-sapply(strsplit(mapped[[2]], ':'), function(x) x[2]);  
    gn<-mapped[,2];
    mp<-split(gn, pth);
    mp<-mp[names(desc)];
    names(mp)<-names(desc);
    mp<-lapply(mp, unique);
    mp<-lapply(mp, function(x) x[!is.na(x)]);
    
    # Save results
    fn<-paste(path, '/r/', nm, '_pathway2gene.rds', sep='');
    saveRDS(list(Organism=species[nm], Pathway=desc, Pathway2Gene=mp), file=fn);
    
    id;
  })
  
  # save log, all pathway IDs
  log.fn<-paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg/log.rds', sep='/');
  if (file.exists(log.fn)) log<-readRDS(log.fn) else log<-list();
  tm<-strsplit(as.character(Sys.time()), ' ')[[1]][1];
  log[[tm]]<-ids;
  saveRDS(log, file=paste((paste(Sys.getenv("RCHIVE_HOME"), 'data/gene.set/public/kegg/log.rds', sep='/'))));
  
  ids;
}
