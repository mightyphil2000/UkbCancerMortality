format_smoking2<-function(dat=NULL){
	dat$smoking<-as.character(dat$smoking)
	dat$smoking[dat$smoking %in% c("Current","Previous")]<-"Ever"	
	dat$smoking[dat$smoking == "Prefer not to answer"]<-NA
	dat$smoking2<-NA
	dat$smoking2[which(dat$smoking=="Ever")]<-1
	dat$smoking2[which(dat$smoking=="Never")]<-0
	return(dat)
}


format_smoking<-function(bd=NULL){
	lvl.0090 <- c(-3,0,1,2)
	lbl.0090 <- c("Prefer not to answer","Never","Previous","Current")
	bd$f.20116.0.0 <- ordered(bd$f.20116.0.0, levels=lvl.0090, labels=lbl.0090)
	bd$f.20116.1.0 <- ordered(bd$f.20116.1.0, levels=lvl.0090, labels=lbl.0090)
	bd$f.20116.2.0 <- ordered(bd$f.20116.2.0, levels=lvl.0090, labels=lbl.0090)
	names(bd)[names(bd) == "f.20116.0.0"]<-"smoking"
	return(bd)
}

add_principal_components<-function(bd=NULL){
  pca<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/principal_components/data.pca1-10.plink.txt",head=FALSE,sep=" ",stringsAsFactors=FALSE)
  names(pca)[names(pca) %in% c("V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")]<-paste0("pc",1:10)
  bd<-merge(bd,pca,by.x=c("geneticID"),by.y="V1")
  return(bd)
}


add_sex_array<-function(bd=NULL){
  covar<-read.table("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/standard_covariates/data.covariates.plink.txt",sep=" ",head=FALSE,stringsAsFactors=FALSE)
  bd<-merge(bd,covar,by.x=c("geneticID"),by.y="V1")
  names(bd)[names(bd) == "V4"]<-"array"
  names(bd)[names(bd) == "V3"]<-"sex2"  
  return(bd)
}
