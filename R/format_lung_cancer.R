#' Format lung cancer mortality data 
#'
#' Format the lung cancer mortality dataset so that it can be used in survival analyses. The 
#' function estimates the number of deaths due to lung cancer and other causes, date of diagnosis, date of death, difference between latest and earlies dates of diagnosis for individuals with more than one diagnosis date, length of followup and survival time.  
#'
#' @param dat the lung cancer dataset containing only individuals with a diagnosis of lung cancer
#' @param censor_date the date up until which the death registry data in UK Biobank is assumed to be complete. We assume that everyone not in the death registry at this date was still alive on this date. The default assumes a censoring date of "2018-02-06", based on the maximum value for the date of diagnosis f.40005.0.0 field in dataset 37205 on the rdsf /projects/MRC-IEU/research/data/ukbiobank/phenotypic/applications/15825/released/2019-05-02/data/derived/formats/r
#'
#' @return data frame
#' @export

format_lung_cancer<-function(dat=NULL,censor_date="2018-02-06"){
	# censor_date in yyyy-mm-dd
	# Secondary cause of death	f.40002
	# Primary cause of death	f.40001
	# f.40000 date of death 
	
	dat<-date_diagnosis(dat=dat) #looks for earliest diagnosis date
	dat<-format_date_of_death(dat=dat)
	dat<-find_max_date_diagnosis(dat=dat)
	dat<-length_followup(dat=dat,censor_date=censor_date)
	dat<-number_deaths_lung_cancer(dat=dat)
	Diff<-as.numeric(dat$max_date_diagnosis - dat$date_diagnosis)
	dat$date_diagnosis_range_days<-Diff
	dat$date_diagnosis_range_months<-round(Diff/(365.25/12),2)
	dat$date_diagnosis_range_years<-round(Diff/(365.25),2)
	return(dat)
}

number_deaths_lung_cancer<-function(dat=NULL){
	Primary<-names(dat)[grep("f.40001",names(dat))]
	Secondary<-names(dat)[grep("f.40002",names(dat))]
	# length(which(!is.na(dat$date_of_death_40000 )))

	Pos0<-which(!is.na(dat$f.40001.0.0 ))
	Pos1<-which(!is.na(dat$f.40001.1.0 ))

	if(!Pos1 %in% Pos0) stop("second instance of f.40001 primary cause of death is not nested within the first instance")
	
	if(!all(dat$f.40001.1.0[Pos1] == dat$f.40001.0.0[Pos1]  )) stop("second instance of f.40001 primary cause of death is different to the first instance")

	dat$primary<-NA
	dat$primary[grep("C",dat$f.40001.0.0)]<-"other malignant neoplasm"
	dat$primary[grep("C34",dat$f.40001.0.0)]<-"Malignant neoplasm of bronchus and lung"
	# dat$primary[grep("C80",dat$f.40001.0.0)]<-"Malignant neoplasm of bronchus and lung" 
	#C80 refers to Malignant neoplasm without specification of site C80-. Probably these are deaths from lung cancer
	Pos<-which(is.na(dat$primary) & !is.na(dat$f.40001.0.0))
	Pos2<-grep("I",dat$f.40001.0.0[Pos])
	dat$primary[Pos[Pos2]] <-"Diseases of the circulatory system"
	Pos<-which(is.na(dat$primary) & !is.na(dat$f.40001.0.0))
	Pos2<-grep("J",dat$f.40001.0.0[Pos])
	dat$primary[Pos[Pos2]] <-"Diseases of the respiratory system"
	Pos<-which(is.na(dat$primary) & !is.na(dat$f.40001.0.0))
	Pos2<-grep("K",dat$f.40001.0.0[Pos])
	dat$primary[Pos[Pos2]] <-"Diseases of the digestive system"
	Pos<-which(is.na(dat$primary) & !is.na(dat$f.40001.0.0))
	dat$primary[Pos] <-"other diseases"

	Names<-names(dat)[grep("f.40001",names(dat))]
	Names2<-names(dat)[grep("f.40002",names(dat))]
	Test<-dat[which( !is.na(dat$date_death) & is.na(dat$primary)),c("date_death","primary","date_diagnosis",Names,Names2)]
	Names<-c(Names,Names2)
	if(!all(is.na(Test[,Names]))) stop("primary cause of death is NA but not all primary and secondary causes of death are NA")
	dat$primary[which( !is.na(dat$date_death) & is.na(dat$primary))]<-"unknown cause of death"
	return(dat)	
}

date_diagnosis<-function(dat=NULL){
		# Pos<-which(!is.na(dat$f.40005.1.0))
		# length(Pos)		
		Names<-names(dat)[grep("f.40005",names(dat))]
		
		x<-1
		N_diag<-lapply(1:nrow(dat),FUN=function(x) 
			unique(as.numeric(dat[x,Names])))

		dat$N_diagnoses<-unlist(lapply(1:nrow(dat),FUN=function(x) 
			length(unlist(N_diag[x])[!is.na(unlist(N_diag[x]))])))

		min_date_diagnosis<-unlist(lapply(1:nrow(dat),FUN=function(x) 
			min(as.numeric(dat[x,Names]),na.rm=TRUE)))
		dat$date_diagnosis<-as.Date(min_date_diagnosis,origin="1970-01-01")

		return(dat)
}
		# bd[Pos[2],c("date_diagnosis",Names)]

format_date_of_death<-function(dat=NULL){
	# Date of death	f.40000 . two instances, the second is nested within, and is identical to, the first 
	Pos0<-which(!is.na(dat$f.40000.0.0))
	Pos1<-which(!is.na(dat$f.40000.1.0))
	if(!all(Pos1 %in% Pos0)) stop("second instance of date of death is not nested within the first instance")
	if(!all(dat$f.40000.1.0[Pos1] == dat$f.40000.0.0[Pos1])) stop("second instance of date of death is different to the first instance")
	dat$f.40000.0.0 <- as.Date(dat$f.40000.0.0)
	dat$f.40000.1.0 <- as.Date(dat$f.40000.1.0)
	names(dat)[names(dat)=="f.40000.0.0"]<-"date_death"
	return(dat)	
}

# time from diagnosis to death or last known date assumed alive
length_followup<-function(dat=NULL,censor_date=NULL){
	dat$survival_days<-as.numeric(dat$date_death -dat$date_diagnosis)
	dat$survival_months<-round(dat$survival_days/(365.25/12),2)
	dat$survival_years<-round(dat$survival_days/365.25,2)

	dat$follow_days<-as.numeric(as.Date(censor_date)-dat$date_diagnosis)
	dat$follow_months<-round(dat$follow_days/(365.25/12),2)
	dat$follow_years<-round(dat$follow_days/365.25,2)
	return(dat)
}


find_max_date_diagnosis<-function(dat=NULL){
	Names<-names(dat)[grep("f.40005",names(dat))]
	max_date_diagnosis<-unlist(lapply(1:nrow(dat),FUN=function(x) 
			max(as.numeric(dat[x,Names]),na.rm=TRUE)))
	dat$max_date_diagnosis<-as.Date(max_date_diagnosis,origin="1970-01-01")
	dat$max_date_diagnosis[which(as.numeric(dat$max_date_diagnosis) == -Inf)]<-NA
	return(dat)
}

