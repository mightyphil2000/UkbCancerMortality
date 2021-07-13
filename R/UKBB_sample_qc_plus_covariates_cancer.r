
restrict_to_white_british<-function(){
  white_british<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/ancestry/data.white_british.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  bd<-merge(bd,white_british,by.x=c("geneticID"),by.y="V1")
  return(bd)
}

add_principal_components<-function(){
  pca<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/principal_components/data.pca1-10.plink.txt",head=FALSE,sep=" ",stringsAsFactors=FALSE)
  names(pca)[names(pca) %in% c("V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")]<-paste0("pc",1:10)
  bd<-merge(bd,pca,by.x=c("geneticID"),by.y="V1")
  return(bd)
}

exclude_relateds<-function(){
  related_high<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/relateds_exclusions/data.highly_relateds.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  related_min<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/relateds_exclusions/data.minimal_relateds.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  bd<-bd[!bd$geneticID %in% related_high$V1, ]
  bd<-bd[!bd$geneticID %in% related_min$V1, ]
  return(bd)
}

standard_exclusion<-function(){
  standard_exclusions<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/standard_exclusions/data.combined_recommended.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  bd<-bd[!bd$geneticID %in% standard_exclusions$V1, ]
  return(bd)
}

add_sex_array<-function(){
  covar<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/standard_covariates/data.covariates.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  bd<-merge(bd,covar,by.x=c("geneticID"),by.y="V1")
  names(bd)[names(bd) == "V4"]<-"array"
  names(bd)[names(bd) == "V3"]<-"sex2"  
  return(bd)
}

# particpant withdrawals
# Dir_extra is for additional lists of widthdrawn participants, e.g. very recenty withdrawals that may not yet have been added to the standard ukb project space on the rdsf 
withdrawn_participants<-function(Dir=NULL,Dir_extra=NULL){
  Files<-paste0(Dir,dir(Dir))
  W<-unlist(lapply(Files,FUN=function(x)
    readLines(x)))

  if(!is.null(Dir_extra)){
    Files<-paste0(Dir_extra,dir(Dir_extra))
    W2<-unlist(lapply(Files,FUN=function(x)
      readLines(x)))
    W<-c(W,W2)
  }
  bd<-bd[!bd$projectID %in% W,]
  return(bd)
}
