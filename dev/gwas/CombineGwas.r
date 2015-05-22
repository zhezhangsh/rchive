# Combine more than one collections of GWAS results
CombineGwas<-function(paths=c()) {
  # paths   Paths to GWAS collections
  
  # names of required files in each collection
  fn.req<-c("analysis2pubmed.rds", "analysis_by_id.rds", "analysis.rds", 
            "pubmed.rds", "study2pubmed.rds", "study_by_id.rds", "study.rds");
}