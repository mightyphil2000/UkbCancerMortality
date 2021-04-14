#' Format cancer mortality data 
#'
#' Format cancer mortality dataset so that it can be used in survival analyses. The 
#' function estimates the number of deaths due to the target cancer of interest and other causes, date of diagnosis, date of death, difference between latest and earlies dates of diagnosis for individuals with more than one diagnosis date, length of followup and survival time. There are two instances of primary cause of death in UK Biobank (f.40001.0.0 and f.40001.1.0). If the second instance is different to the first instance (ie an individual has more than one primary cause of death), the function checks to see if either of these instances match the indicated cancer. If there is a match, then the primary cause of death is set to that cancer. E.g. one a case has prostate cancer (C61) and Pulmonary edema (J81) listed as primary causes of death. In this case, we would set the primary cause of death to prostate cancer. This is a rare occurence. I've only observed this example once for prostate cancer and 0 times for lung cancer, breast cancer and overall cancer. If a case has two primary causes of death and neither one matches the indicated cancer, then the function stops and returns a warning message. 
#'
#' @param dat the cancer dataset 
#' @param censor_date the date up until which the death registry data in UK Biobank is assumed to be complete. We assume that everyone not in the death registry at this date was still alive on this date. The default assumes a censoring date of "2018-02-06", based on the maximum value for the date of diagnosis f.40005.0.0 field in dataset 37205 on the rdsf /projects/MRC-IEU/research/data/ukbiobank/phenotypic/applications/15825/released/2019-05-02/data/derived/formats/r
#' @param cancer_name the name of the field/column in dat corresponding to case-control status for the cancer of interest. The field must take on values of 1 (control) or 2 (cases). 
#' @param icd10 the icd10 code for the cancer of interest. Icd codes are hierarchically structured. Therefore, it is not necessary to include all the lowest level icd codes for the cancer of interest. Only the highest level code should be specified. For example, breast cancer corresponds to icd10 codes C50.0-C50.9. Therefore, to create a breast cancer dataset, you would set icd10 to "C50". The function then selects all icd10 codes containing the string "C50". 
#'
#' @return data frame
#' @export

format_cancer<-function(dat=NULL,censor_date="2018-02-06",cancer_name=NULL,icd10=NULL){
	# censor_date in yyyy-mm-dd
	# Secondary cause of death	f.40002
	# Primary cause of death	f.40001
	# f.40000 date of death 

	dat<-dat[which(dat[,cancer_name] == 2),]	
	dat<-date_diagnosis(dat=dat) #looks for earliest diagnosis date
	dat<-format_date_of_death(dat=dat)
	dat<-find_max_date_diagnosis(dat=dat)	
	dat<-length_followup(dat=dat,censor_date=censor_date)
	# dat2<-dat
	dat<-number_deaths_cancer(dat=dat,icd10=icd10)
	Diff<-as.numeric(dat$max_date_diagnosis - dat$date_diagnosis)
	dat$date_diagnosis_range_days<-Diff
	dat$date_diagnosis_range_months<-round(Diff/(365.25/12),2)
	dat$date_diagnosis_range_years<-round(Diff/(365.25),2)
	return(dat)
}

number_deaths_cancer<-function(dat=NULL,icd10=icd10){
	Primary<-names(dat)[grep("f.40001",names(dat))]
	Secondary<-names(dat)[grep("f.40002",names(dat))]
	# length(which(!is.na(dat$date_of_death_40000 )))

	Pos0<-which(!is.na(dat$f.40001.0.0 ))
	Pos1<-which(!is.na(dat$f.40001.1.0 ))
	
	# dat$f.40001.0.0[Pos1]
	# dat$f.40001.1.0[Pos1]
	# dat$f.40002.0.0[Pos1]
	# dat$f.40002.1.0[Pos1]
	# dat$f.40002.2.0[Pos1]

	if(!all(Pos1 %in% Pos0)) stop("second instance of f.40001 primary cause of death is not nested within the first instance")
	
	cancer_name2<-gsub("_"," ",cancer_name)
	cancer_name2<-gsub("overall","",cancer_name2)
	cancer_name2<-gsub("cancer","",cancer_name2)
	cancer_name2<-trimws(cancer_name2)
	cancer_name2<-paste0("Malignant neoplasm of ",cancer_name2)
	dat$primary<-rep(NA,nrow(dat))

	if(!all(dat$f.40001.1.0[Pos1] == dat$f.40001.0.0[Pos1]  )) {
		# position where first instance does not equal second instance of f.40001
		Pos1_2<-Pos1[which(dat$f.40001.1.0[Pos1] != dat$f.40001.0.0[Pos1] )]
		Dat2<-dat[Pos1_2,]
		Dat1<-dat[1:nrow(dat) != Pos1_2,]
		Pos1_3<-which(Dat2$f.40001.1.0 ==icd10 | Dat2$f.40001.0.0 ==icd10)
		# if either the first or second instance match the cancer of interest, set primary cause of death to the cancer of interest. This means that when a case has two different primary causes of death, one of which matches the cancer of interest, we set primary cause of death to the cancer of interest. E.g. one case had prostate cancer (C61) and Pulmonary edema (J81) listed as primary causes of death. We set the primary cause of death for this case to prostate cancer. This is a rare occurence. I've only observed this example once for prostate cancer and 0 times for lung cancer, breast cancer and overall cancer
		Dat2$primary[Pos1_3]<-cancer_name2		
		if(Dat2$f.40001.1.0 !=icd10 & Dat2$f.40001.0.0 !=icd10) stop("second instance of f.40001 primary cause of death is different to the first instance and neither corresponds to the indicated cancer of interest")
		dat<-rbind(Dat1,Dat2)
	}	

	
	dat$primary[grep("C",dat$f.40001.0.0)]<-"other malignant neoplasm"
	dat$primary[grep(icd10,dat$f.40001.0.0)]<-cancer_name2
	
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
	Pos2<-grep("G",dat$f.40001.0.0[Pos])
	dat$primary[Pos[Pos2]] <-"Diseases of the nervous system"
	Pos<-which(is.na(dat$primary) & !is.na(dat$f.40001.0.0))
	Pos2<-unique(unlist(lapply(c("D0","D1","D2","D3","D4"),FUN=function(x)
		grep(x,dat$f.40001.0.0[Pos]))))
	dat$primary[Pos[Pos2]] <-"In situ or benign neoplasms or neoplasms of unknown/uncertain behavior, polycythemia vera and myelodysplastic syndromes"
	Pos<-which(is.na(dat$primary) & !is.na(dat$f.40001.0.0))
	Pos2<-grep("F",dat$f.40001.0.0[Pos])
	dat$primary[Pos[Pos2]] <-"Mental, Behavioral and Neurodevelopmental disorders"
	Pos<-which(is.na(dat$primary) & !is.na(dat$f.40001.0.0))
	Pos2<-unique(unlist(lapply(c("A","B"),FUN=function(x)
		grep(x,dat$f.40001.0.0[Pos]))))
	dat$primary[Pos[Pos2]] <-"Certain infectious and parasitic diseases"
	
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