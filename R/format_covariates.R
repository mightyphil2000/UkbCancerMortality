format_smoking<-function(dat=NULL){
	dat$smoking<-as.character(dat$smoking)
	dat$smoking[dat$smoking %in% c("Current","Previous")]<-"Ever"	
	dat$smoking[dat$smoking == "Prefer not to answer"]<-NA
	dat$smoking2<-NA
	dat$smoking2[which(dat$smoking=="Ever")]<-1
	dat$smoking2[which(dat$smoking=="Never")]<-0
	return(dat)
}

