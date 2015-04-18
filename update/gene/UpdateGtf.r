library(devtools);
source_url("https://raw.githubusercontent.com/zhezhangsh/rchive/master/load.r");

path=paste(RCHIVE_HOME, 'data/gene/public/gtf', sep='/');

if (!file.exists(path)) dir.create(path, recursive=TRUE);
if(!file.exists(paste(path, 'r', sep='/'))) dir.create(paste(path, 'r', sep='/'), recursive=TRUE);
if(!file.exists(paste(path, 'src', sep='/'))) dir.create(paste(path, 'src', sep='/'), recursive=TRUE);

##############################################################################################################
# Files to be downloaed directly from sources
fn<-c(
  #'human_GRCh38_NCBI_Scaffold' = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p2_scaffolds.gff3.gz", 
  'human_GRCh38_NCBI' = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p2_top_level.gff3.gz",
  #'human_GRCh37_NCBI_Scaffold' = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.3/GFF/ref_GRCh37.p5_scaffolds.gff3.gz",
  'human_GRCh37_NCBI' = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.3/GFF/ref_GRCh37.p5_top_level.gff3.gz",
  
  'human_GRCh38_Gencode' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz",
  #'human_GRCh38_Gencode_Scaffold' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_22/gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz",
  'human_GRCh38_Gencode_PolyA' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_22/gencode.v22.polyAs.gtf.gz",
  'human_GRCh38_Gencode_Pseudo' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_22/gencode.v22.2wayconspseudos.gtf.gz",
  'human_GRCh38_Gencode_tRNA' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_22/gencode.v22.tRNAs.gtf.gz",
  'human_GRCh37_Gencode' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
  #'human_GRCh37_Gencode_Scaffold' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz",
  'human_GRCh37_Gencode_PolyA' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.polyAs.gtf.gz",
  'human_GRCh37_Gencode_Pseudo' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.2wayconspseudos.gtf.gz",
  'human_GRCh37_Gencode_tRNA' = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.tRNAs.gtf.gz",

  'human_GRCh38_Ensembl' = "ftp://ftp.ensembl.org/pub/release-79/gtf/homo_sapiens//Homo_sapiens.GRCh38.79.gtf.gz",
  'human_GRCh37_Ensembl' = "ftp://ftp.ensembl.org/pub//release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz",
  
  'human_GRCh37_AceView' = "ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_37_Aug10.human.genes/AceView.ncbi_37.genes_gff.gff.gz",
  
  
  #'mouse_GRCm38_NCBI_Scaffold' = "ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/GFF/ref_GRCm38.p3_scaffolds.gff3.gz",
  'mouse_GRCm38_NCBI' = "ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/GFF/ref_GRCm38.p3_top_level.gff3.gz",
  #'mouse_MGSCv37_NCBI_Scaffold' = "ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/ARCHIVE/BUILD.37.2/GFF/ref_MGSCv37_scaffolds.gff3.gz",
  'mouse_MGSCv37_NCBI' = "ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/ARCHIVE/BUILD.37.2/GFF/ref_MGSCv37_top_level.gff3.gz",
  'mouse_MGSCv37_AceView' = "ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_37_Sep07.mouse.genes/AceView.mm_37.genes_gff.tar.gz",
  'mouse_GRCm38_Ensembl' = "ftp://ftp.ensembl.org/pub/release-79/gtf//mus_musculus/Mus_musculus.GRCm38.79.gtf.gz",

  'worm_WB235_Ensembl' = "ftp://ftp.ensembl.org/pub/release-79/gtf//caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.79.gtf.gz"
)

fn0<-sapply(strsplit(fn, '/'), function(x) x[length(x)]);
fn0<-paste(path, 'src', fn0, sep='/');
names(fn0)<-names(fn);
fn1<-sapply(names(fn), function(nm) if(!file.exists(fn0[nm])) download.file(fn[nm], fn0[nm]));

##############################################################################################################
# GTF files previously downloaded from UCSC using Table Browser
fn2<-dir(paste(path, 'src', sep='/'));
fn2<-fn2[grep('UCSC', fn2)];
nm<-sapply(strsplit(fn2, '\\.'), function(x) x[1]);
fn2<-paste(path, 'src', fn2, sep='/');
names(fn2)<-nm;

fn.gtf<-c(fn0, fn2); # All source gtf files

###################
fn<-sapply(names(fn.gtf), function(nm) {
  fn<-fn.gtf[]
})
