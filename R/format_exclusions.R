# Exclude individuals with multiple cancer diagnoses if target cancer (e.g. LC) is first cancer diagnosis but a non-target cancer diagnosis occurs within 6 months of the first diagnosis. The concern here is that the target cancer diagnosis is not the primary cancer but has spread to the target site from a different primary site
exclude_6mo<-function(bd=NULL){
	bd<-bd[!bd$Index_6mo,]
	return(bd)
}

# Exclude individuals with with multiple cancer diagnoses if target cancer (eg LC) is not first diagnosis and occurs <3years after the first diagnosis. The concern here is that the target cancer is actually a relapse of a previously diagnosed non-target cancer
exclude_3yr<-function(bd=NULL){
  bd<-bd[!bd$Index_3yr,]
  return(bd)
}

# 5.Exclude individuals who die within 3 months of a diagnosis 
# the concern here is that cases that die within 3 months are late stage diagnoses, which may reduce the heritability of lung cancer survival 
exclude_death_3months<-function(bd=NULL){
  bd<-bd[which(is.na(bd$survival_months) | bd$survival_months>=3),]
  return(bd)
}


# particpant withdrawals
# Dir_extra is for additional lists of widthdrawn participants, e.g. very recenty withdrawals that may not yet have been added to the standard ukb project space on the rdsf 
withdrawn_participants<-function(bd=NULL,Dir=NULL,Dir_extra=NULL){
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


exclude_relateds<-function(bd=NULL){
  related_high<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/relateds_exclusions/data.highly_relateds.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  related_min<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/relateds_exclusions/data.minimal_relateds.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  bd<-bd[!bd$geneticID %in% related_high$V1, ]
  bd<-bd[!bd$geneticID %in% related_min$V1, ]
  return(bd)
}

standard_exclusion<-function(bd=NULL){
  standard_exclusions<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/standard_exclusions/data.combined_recommended.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  bd<-bd[!bd$geneticID %in% standard_exclusions$V1, ]
  return(bd)
}




restrict_to_white_british<-function(bd=NULL){
  white_british<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/ancestry/data.white_british.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  bd<-merge(bd,white_british,by.x=c("geneticID"),by.y="V1")
  return(bd)
}

# 6.exclude prevalent cases - 312 prevalent cases, 2517 incident cases . Left truncation worth the hassle? - done
exclude_prevalent_cases<-function(bd=NULL){
  bd<-bd[which(bd$incident_lung_cancer==2),]
  return(bd)
}
