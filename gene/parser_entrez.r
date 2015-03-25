# Parse the NCBI Entrez genes of multiple species
path<-paste(RCHIVE_HOME, 'gene/public/entrez', sep='/');

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'));
if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'));

# Download Entrez genes from NCBI FTP server
fn.src<-c(
    'human' = 'ftp://ftp.ncbi.nih.gov/gene//DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz', 
    'mouse' = 'ftp://ftp.ncbi.nih.gov/gene//DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz',
    'rat' = 'ftp://ftp.ncbi.nih.gov/gene//DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz',
    'zebrafish' = 'ftp://ftp.ncbi.nih.gov/gene//DATA/GENE_INFO/Non-mammalian_vertebrates/Danio_rerio.gene_info.gz',
    'fly' = 'ftp://ftp.ncbi.nih.gov/gene//DATA/GENE_INFO/Invertebrates/Drosophila_melanogaster.gene_info.gz',
    'worm' = 'ftp://ftp.ncbi.nih.gov/gene//DATA/GENE_INFO/Invertebrates/Caenorhabditis_elegans.gene_info.gz'
)

    ############################################################
    # Parse the NCBI Entrez genes of a species
    parser<-function(fn.src, species, path) {
        if (!require(plyr)) {
          install.packages('plyr');
          require(plyr);
        };
                
        fn<-paste(path, 'src', species, sep='/');
        download.file(fn.src, fn); # download file from NCBI
        
        # read in all info.
        raw<-scan(fn, sep='\n', what='', flush=TRUE);
        
        # header line
        hd<-strsplit(sub('^#Format: ', '', raw[1]), ' ')[[1]];
        
        # create final annotation table
        gn<-data.frame(do.call('rbind', strsplit(raw[-1], '\t')));
        names(gn)<-hd[1:ncol(gn)];
        rownames(gn)<-as.vector(gn$GeneID); # use Gene ID as row names
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
        saveRDS(gn, file=paste(path, '/r/', species, '_genes_full.rds', sep=''));
        saveRDS(gn.tbl, file=paste(path, '/r/', species, '_genes_datatable.rds', sep=''));
        saveRDS(id2sym, file=paste(path, '/r/', species, '_genes_id2symbol.rds', sep=''));
        saveRDS(sym2id, file=paste(path, '/r/', species, '_genes_symbol2id.rds', sep=''));
        saveRDS(id2syn, file=paste(path, '/r/', species, '_genes_id2synonyms.rds', sep=''));
        saveRDS(syn2id, file=paste(path, '/r/', species, '_genes_synonyms2id.rds', sep=''));
    }
    ############################################################

# All genes IDs after update
fns<-paste(path, '/r/', names(fn.src), '_genes_full.rds', sep='');
ids.old<-lapply(fns, function(f) if (file.exists(f)) rownames(readRDS(f)) else c());
names(ids.old)<-names(fn.src);

for (i in 1:length(fn.src)) {
    print(names(fn.src)[i]);
    parser(fn.src[i], names(fn.src)[i], path);
}

# All genes IDs after update
ids.new<-lapply(fns, function(f) rownames(readRDS(f)));
names(ids.new)<-names(fn.src);
up<-list(N=sapply(ids.new, length), 
         Added=lapply(names(ids.new), function(i) setdiff(ids.new[[i]], ids.old[[i]])), 
         Removed=lapply(names(ids.old), function(nm) setdiff(ids.old[[i]], ids.new[[i]])));

log<-readRDS(paste(path, 'log.rds', sep='/'));
log<-c(log, list(up));
names(log)[length(log)]<-as.character(Sys.Date());
saveRDS(log, file=paste(path, 'log.rds', sep='/'));