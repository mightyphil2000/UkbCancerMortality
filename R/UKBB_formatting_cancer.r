#' @importFrom magrittr %>%

lung_cancer_function<-function(){
	#lung cancer
	require(tidyr)
	require(dplyr)
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C340", "C341", "C342", "C343", "C348", "C349")
	lungICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("1622", "1623", "1624", "1625", "1628", "1629")
	lungICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	lungICD9 <-as.integer(lungICD9)
	

	bd<-UKBcancerFunc(dat=bd,cancerCode = lungICD9,sitename = "lungICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = lungICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = lungICD10,sitename = "lungICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = lungICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent
	bd$lungPrevelent<-ifelse(!is.na(bd$lungICD101), 1,NA)
	bd$lungPrevelent<-ifelse(!is.na(bd$lungICD91), 1,bd$lungPrevelent)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)
	bd$lungSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
							ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
												ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	
	return(bd)
}



acute_lymphoblastic_leukemia_function<-function(icd9=NULL,icd10=NULL){
	# acute lymphoblastic leukemia 
		
	#acute_lymph_leuk cancer
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C910")
	acute_lymph_leukICD10 <- unique (grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(acute_lymph_leukICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("2040")
	acute_lymph_leukICD9 <- unique (grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	acute_lymph_leukICD9 <-as.integer(acute_lymph_leukICD9)
	table(acute_lymph_leukICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = acute_lymph_leukICD9,sitename = "acute_lymph_leukICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = acute_lymph_leukICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = acute_lymph_leukICD10,sitename = "acute_lymph_leukICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = acute_lymph_leukICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$acute_lymph_leukPrevelent<-ifelse(!is.na(bd$acute_lymph_leukICD101), 1,NA)
	bd$acute_lymph_leukPrevelent<-ifelse(!is.na(bd$acute_lymph_leukICD91), 1,bd$acute_lymph_leukPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$acute_lymph_leukSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	return(bd)
}

endometrial_cancer_function<-function(icd9=NULL,icd10=NULL){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)

	# toMatch <- c("C530","C531","C538","C539")
	# allICD10[grep("C541",allICD10)]
	# toMatch <-  c("C641","C642","C649") 
	toMatch <-  c( "C541" ) 
	endometrialICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(endometrialICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]


	toMatch <- c(1820 )
	endometrialICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	endometrialICD9 <-as.integer(endometrialICD9)
	table(endometrialICD9)

	# allICD9[grep(2050,allICD9)]

	bd<-UKBcancerFunc(dat=bd,cancerCode = endometrialICD9,sitename = "endometrialICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = endometrialICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = endometrialICD10,sitename = "endometrialICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = endometrialICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent
	bd$endometrialPrevelent<-ifelse(!is.na(bd$endometrialICD101), 1,NA)
	bd$endometrialPrevelent<-ifelse(!is.na(bd$endometrialICD91), 1,bd$endometrialPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report endometrials...")

	#self report
	#identify self reported all endometrial (coded as 1)
	bd$endometrialSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	 # rm(dat)
	return(bd)
}

# bd4<-bd
aml_cancer_function<-function(icd9=NULL,icd10=NULL){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)

	# toMatch <- c("C530","C531","C538","C539")
	# allICD10[grep("C920",allICD10)]
	# toMatch <-  c("C641","C642","C649") 
	toMatch <-  c( "C920" ) 
	amlICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(amlICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]


	toMatch <- c(2050 )
	amlICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	amlICD9 <-as.integer(amlICD9)
	table(amlICD9)

	# allICD9[grep(2050,allICD9)]

	bd<-UKBcancerFunc(dat=bd,cancerCode = amlICD9,sitename = "amlICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = amlICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = amlICD10,sitename = "amlICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = amlICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent
	bd$amlPrevelent<-ifelse(!is.na(bd$amlICD101), 1,NA)
	bd$amlPrevelent<-ifelse(!is.na(bd$amlICD91), 1,bd$amlPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report amls...")

	#self report
	#identify self reported all aml (coded as 1)
	bd$amlSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	 # rm(dat)
	return(bd)
}

# bd4<-bd
stomach_cancer_function<-function(icd9=NULL,icd10=NULL){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)

	# toMatch <- c("C530","C531","C538","C539")
	# allICD10[grep("C16",allICD10)]
	# toMatch <-  c("C641","C642","C649") 
	toMatch <-  c( "C160","C161" ,"C162","C163","C164","C165",  "C166", "C168","C169" ) 
	stomachICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(stomachICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]


	toMatch <- c(1510, 1512, 1514,   1515,1519 )
	stomachICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	stomachICD9 <-as.integer(stomachICD9)
	table(stomachICD9)

	allICD9[grep(151,allICD9)]

	bd<-UKBcancerFunc(dat=bd,cancerCode = stomachICD9,sitename = "stomachICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = stomachICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = stomachICD10,sitename = "stomachICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = stomachICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent
	bd$stomachPrevelent<-ifelse(!is.na(bd$stomachICD101), 1,NA)
	bd$stomachPrevelent<-ifelse(!is.na(bd$stomachICD91), 1,bd$stomachPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report stomachs...")

	#self report
	#identify self reported all stomach (coded as 1)
	bd$stomachSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	 # rm(dat)
	return(bd)
}

# bd4<-bd
pancreatic_cancer_function<-function(icd9=NULL,icd10=NULL){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)

	# toMatch <- c("C530","C531","C538","C539")
	# allICD10[grep("C25",allICD10)]
	# toMatch <-  c("C641","C642","C649") 
	toMatch <-  c( "C250", "C251", "C252", "C253","C254","C257" ,"C258" ,"C259") 
	pancreaticICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(pancreaticICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]


	toMatch <- c(1570,1571,1572,1573,1574,1578,1579)
	pancreaticICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	pancreaticICD9 <-as.integer(pancreaticICD9)
	table(pancreaticICD9)

	# allICD9[grep(157,allICD9)]

	bd<-UKBcancerFunc(dat=bd,cancerCode = pancreaticICD9,sitename = "pancreaticICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = pancreaticICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = pancreaticICD10,sitename = "pancreaticICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = pancreaticICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent
	bd$pancreaticPrevelent<-ifelse(!is.na(bd$pancreaticICD101), 1,NA)
	bd$pancreaticPrevelent<-ifelse(!is.na(bd$pancreaticICD91), 1,bd$pancreaticPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report pancreatics...")

	#self report
	#identify self reported all pancreatic (coded as 1)
	bd$pancreaticSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	 # rm(dat)
	return(bd)
}



# bd1<-bd
kidney_cancer_function<-function(icd9=NULL,icd10=NULL){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)

	# toMatch <- c("C530","C531","C538","C539")
	# allICD10[grep("64",allICD10)]
	# toMatch <-  c("C641","C642","C649") 
	toMatch <-  c("C64") 
	kidneyICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(kidneyICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- 1890 
	kidneyICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	kidneyICD9 <-as.integer(kidneyICD9)
	table(kidneyICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = kidneyICD9,sitename = "kidneyICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = kidneyICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = kidneyICD10,sitename = "kidneyICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = kidneyICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent
	bd$kidneyPrevelent<-ifelse(!is.na(bd$kidneyICD101), 1,NA)
	bd$kidneyPrevelent<-ifelse(!is.na(bd$kidneyICD91), 1,bd$kidneyPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report kidneys...")

	#self report
	#identify self reported all kidney (coded as 1)
	bd$kidneySelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	 # rm(dat)
	return(bd)
}


cervical_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)

	toMatch <- c("C530","C531","C538","C539")
	cervicalICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(cervicalICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c(1800,1801,1808,1809)
	cervicalICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	cervicalICD9 <-as.integer(cervicalICD9)
	table(cervicalICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = cervicalICD9,sitename = "cervicalICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = cervicalICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = cervicalICD10,sitename = "cervicalICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = cervicalICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent
	bd$cervicalPrevelent<-ifelse(!is.na(bd$cervicalICD101), 1,NA)
	bd$cervicalPrevelent<-ifelse(!is.na(bd$cervicalICD91), 1,bd$cervicalPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$cervicalSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	 # rm(dat)
	return(bd)
}

overall_cancer_function<-function(){
	#All cancer
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)

	# take out skin cancer
	#allICD10 <- grep("C44", allICD10, value=T, invert=T)


	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	x <-c("140", "141", "142", "143", "144", "145", "146", "147", "148", "149", "150", "151", "152", "153", "154", "155", "156", "157", "158", "159", "160", "161", "162", "163", "164", "165", "166", "167", "168", "169", "170", "171", "172", "173", "174", "175", "176", "177", "178", "179", "180", "181", "182", "183", "184", "185", "186", "187", "188", "189", "190", "191", "192", "193", "194", "195", "196", "197", "198", "199", "200", "201", "202", "203", "204", "205", "206", "207", "208")
	allICD9 <- grep(paste0(x, collapse = "|"), as.character(allICD9), value = TRUE)
	allICD9 <-as.integer(allICD9)
	# take out skin cancer
	#allICD9 <- grep("173", allICD9, value=T, invert=T)
	
		# unique(bd$f.40012.0.0)
	# unique(bd$f.40012.0.0)

	names(bd)
	
	bd<-UKBcancerFunc(dat=bd,cancerCode = allICD9,sitename = "allICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = allICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = allICD10,sitename = "allICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = allICD10,sitename = "otherICD10", other=T)
	# define prevalent
	bd$allPrevelent<-ifelse(!is.na(bd$allICD101), 1,NA)
	bd$allPrevelent<-ifelse(!is.na(bd$allICD91), 1,bd$allPrevelent)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)
	bd$allSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))

	 # rm(dat)
# df <- as.data.frame(df)
	# df<-bd
	return(bd)
}

# exclude non-melanoma skin cancer
overall_cancer_exclc44_function<-function(){
	# get all ICD10 cancer codes
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls as these are carcinoma in situs)
	allICD10 <- grep("C", allICD10, value=T)
	# take out non-melanoma skin cancer code (C44)
	allICD10 <- grep("C44", allICD10, value=T, invert=T)
	
	# get all ICD9 cancer codes
	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	# keep all relevent ICD9 cancer codes (do not keep carcinoma in situ codes)
	x <-c("140", "141", "142", "143", "144", "145", "146", "147", "148", "149", "150", "151", "152", "153", "154", "155", "156", "157", "158", "159", "160", "161", "162", "163", "164", "165", "166", "167", "168", "169", "170", "171", "172", "173", "174", "175", "176", "177", "178", "179", "180", "181", "182", "183", "184", "185", "186", "187", "188", "189", "190", "191", "192", "193", "194", "195", "196", "197", "198", "199", "200", "201", "202", "203", "204", "205", "206", "207", "208")
	allICD9 <- grep(paste0(x, collapse = "|"), as.character(allICD9), value = TRUE)
	allICD9 <-as.integer(allICD9)
	# take out non-melanoma skin cancer code (173)
	allICD9 <- grep("173", allICD9, value=T, invert=T)

	# Execute UKBcancerFunc function to extract all the cancer phenotypes of interest
	bd<-UKBcancerFunc(dat=bd,cancerCode = allICD9,sitename = "allICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = allICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = allICD10,sitename = "allICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = allICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent - this is actually just creating a variable that says if an individual has cancer or not (at any timepoint)
	bd$allPrevelent<-ifelse(!is.na(bd$allICD101), 1,NA)
	bd$allPrevelent<-ifelse(!is.na(bd$allICD91), 1,bd$allPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	# self report
	# identify self reported all cancers (coded as 1) - these will be used to exclude participants from controls
	bd$allSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)	
}

	# df_split<-separate_cancers()
	# df_split<-format_columns()
	# df_split<-generate_incident_flag()
	# df_split<-generate_behaviour_flag()
	# df_split<-generate_controls()
	# df_split<-generate_incident_cases()
	# bc<-tidy_up()
brain_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C71")
	brainICD10 <- unique (grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(brainICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("191")
	brainICD9 <- unique (grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	brainICD9 <-as.integer(brainICD9)
	table(brainICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = brainICD9,sitename = "brainICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = brainICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = brainICD10,sitename = "brainICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = brainICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$brainPrevelent<-ifelse(!is.na(bd$brainICD101), 1,NA)
	bd$brainPrevelent<-ifelse(!is.na(bd$brainICD91), 1,bd$brainPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$brainSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

breast_cancer_function<-function(){

		allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C500", "C501", "C502", "C503", "C504", "C505", "C506", "C508", "C509")
	breastICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(breastICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c(1740,1741,1742,1743,1744,1745,1746,1747,1748,1749)
	breastICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	breastICD9 <-as.integer(breastICD9)
	table(breastICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = breastICD9,sitename = "breastICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = breastICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = breastICD10,sitename = "breastICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = breastICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define prevalent
	bd$breastPrevelent<-ifelse(!is.na(bd$breastICD101), 1,NA)
	bd$breastPrevelent<-ifelse(!is.na(bd$breastICD91), 1,bd$breastPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$breastSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}
	
melanoma_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C430", "C431", "C432", "C433", "C434", "C435", "C436", "C437", "C438", "C439")
	skinICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(skinICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c(1720, 1721, 1722, 1723, 1724, 1725, 1726, 1727, 1728, 1729)
	skinICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	skinICD9 <-as.integer(skinICD9)
	table(skinICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = skinICD9,sitename = "skinICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = skinICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = skinICD10,sitename = "skinICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = skinICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$skinPrevelent<-ifelse(!is.na(bd$skinICD101), 1,NA)
	bd$skinPrevelent<-ifelse(!is.na(bd$skinICD91), 1,bd$skinPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$skinSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}



prostate_cancer_function<-function(){

	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)

	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C61")
	prostateICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(prostateICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c(185)
	prostateICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	prostateICD9 <-as.integer(prostateICD9)
	table(prostateICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = prostateICD9,sitename = "prostateICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = prostateICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = prostateICD10,sitename = "prostateICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = prostateICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$prostatePrevelent<-ifelse(!is.na(bd$prostateICD101), 1,NA)
	bd$prostatePrevelent<-ifelse(!is.na(bd$prostateICD91), 1,bd$prostatePrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$prostateSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

pharynx_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C01", "C024", "C051", "C052", "C058", "C059", "C090", "C091", "C098", "C099", "C100", "C101", "C102", "C103", "C104", "C108", "C109", "C12", "C130", "C131", "C132", "C139", "C140", "C142")
	pharynxICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(pharynxICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("1410", "1453", "1455", "1460", "1461")
	pharynxICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	pharynxICD9 <-as.integer(pharynxICD9)
	table(pharynxICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = pharynxICD9,sitename = "pharynxICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = pharynxICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = pharynxICD10,sitename = "pharynxICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = pharynxICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$pharynxPrevelent<-ifelse(!is.na(bd$pharynxICD101), 1,NA)
	bd$pharynxPrevelent<-ifelse(!is.na(bd$pharynxICD91), 1,bd$pharynxPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$pharynxSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

ovarian_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C56")
	ovarianICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(ovarianICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c(1830)
	ovarianICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	ovarianICD9 <-as.integer(ovarianICD9)
	table(ovarianICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = ovarianICD9,sitename = "ovarianICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = ovarianICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = ovarianICD10,sitename = "ovarianICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = ovarianICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$ovarianPrevelent<-ifelse(!is.na(bd$ovarianICD101), 1,NA)
	bd$ovarianPrevelent<-ifelse(!is.na(bd$ovarianICD91), 1,bd$ovarianPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$ovarianSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

oropharyngeal_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C003", "C004", "C005", "C006", "C009", "C020", "C021", "C022", "C023", "C028", "C029", "C030", "C031", "C039", "C040", "C041", "C048", "C049", "C050", "C060", "C061", "C062", "C068", "C069", "C01", "C024", "C051", "C052", "C058", "C059", "C090", "C091", "C098", "C099", "C100", "C101", "C102", "C103", "C104", "C108", "C109", "C12", "C130", "C131", "C132", "C139", "C140", "C142")
	oral_pharynxICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(oral_pharynxICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("1412", "1413", "1419", "1430", "1431", "1449", "1450", "1451", "1452")
	oral_pharynxICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	oral_pharynxICD9 <-as.integer(oral_pharynxICD9)
	table(oral_pharynxICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = oral_pharynxICD9,sitename = "oral_pharynxICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = oral_pharynxICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = oral_pharynxICD10,sitename = "oral_pharynxICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = oral_pharynxICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$oral_pharynxPrevelent<-ifelse(!is.na(bd$oral_pharynxICD101), 1,NA)
	bd$oral_pharynxPrevelent<-ifelse(!is.na(bd$oral_pharynxICD91), 1,bd$oral_pharynxPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$oral_pharynxSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

oral_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C003", "C004", "C005", "C006", "C009", "C020", "C021", "C022", "C023", "C028", "C029", "C030", "C031", "C039", "C040", "C041", "C048", "C049", "C050", "C060", "C061", "C062", "C068", "C069")
	oral_cavityICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(oral_cavityICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("1412", "1413", "1419", "1430", "1431", "1449", "1450", "1451", "1452")
	oral_cavityICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	oral_cavityICD9 <-as.integer(oral_cavityICD9)
	table(oral_cavityICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = oral_cavityICD9,sitename = "oral_cavityICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = oral_cavityICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = oral_cavityICD10,sitename = "oral_cavityICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = oral_cavityICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$oral_cavityPrevelent<-ifelse(!is.na(bd$oral_cavityICD101), 1,NA)
	bd$oral_cavityPrevelent<-ifelse(!is.na(bd$oral_cavityICD91), 1,bd$oral_cavityPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$oral_cavitySelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

oesophageal_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C15")
	oesophICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(oesophICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("150")
	oesophICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	oesophICD9 <-as.integer(oesophICD9)
	table(oesophICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = oesophICD9,sitename = "oesophICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = oesophICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = oesophICD10,sitename = "oesophICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = oesophICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$oesophPrevelent<-ifelse(!is.na(bd$oesophICD101), 1,NA)
	bd$oesophPrevelent<-ifelse(!is.na(bd$oesophICD91), 1,bd$oesophPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$oesophSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

melanoma_plus_other_malignant_skin_cancer_function<-function(){
	# Melanoma and other malignant neoplasms of skin
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C43","C44")
	mmplus_skinICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(mmplus_skinICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("172","173")
	mmplus_skinICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	mmplus_skinICD9 <-as.integer(mmplus_skinICD9)
	table(mmplus_skinICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = mmplus_skinICD9,sitename = "mmplus_skinICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = mmplus_skinICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = mmplus_skinICD10,sitename = "mmplus_skinICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = mmplus_skinICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$mmplus_skinPrevelent<-ifelse(!is.na(bd$mmplus_skinICD101), 1,NA)
	bd$mmplus_skinPrevelent<-ifelse(!is.na(bd$mmplus_skinICD91), 1,bd$mmplus_skinPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$mmplus_skinSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}


nonmelanoma_skin_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C44")
	nm_skinICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(nm_skinICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("173")
	nm_skinICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	nm_skinICD9 <-as.integer(nm_skinICD9)
	table(nm_skinICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = nm_skinICD9,sitename = "nm_skinICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = nm_skinICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = nm_skinICD10,sitename = "nm_skinICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = nm_skinICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$nm_skinPrevelent<-ifelse(!is.na(bd$nm_skinICD101), 1,NA)
	bd$nm_skinPrevelent<-ifelse(!is.na(bd$nm_skinICD91), 1,bd$nm_skinPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$nm_skinSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

bladder_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C670", "C671", "C672", "C673", "C674", "C675", "C676", "C677", "C678", "C679")
	bladderICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(bladderICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c(1880, 1882, 1884, 1886, 1888, 1889)
	bladderICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	bladderICD9 <-as.integer(bladderICD9)
	table(bladderICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = bladderICD9,sitename = "bladderICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = bladderICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = bladderICD10,sitename = "bladderICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = bladderICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$bladderPrevelent<-ifelse(!is.na(bd$bladderICD101), 1,NA)
	bd$bladderPrevelent<-ifelse(!is.na(bd$bladderICD91), 1,bd$bladderPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$bladderSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

myeloid_leukemia_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C92")
	myel_leukICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(myel_leukICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("205")
	myel_leukICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	myel_leukICD9 <-as.integer(myel_leukICD9)
	table(myel_leukICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = myel_leukICD9,sitename = "myel_leukICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = myel_leukICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = myel_leukICD10,sitename = "myel_leukICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = myel_leukICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$myel_leukPrevelent<-ifelse(!is.na(bd$myel_leukICD101), 1,NA)
	bd$myel_leukPrevelent<-ifelse(!is.na(bd$myel_leukICD91), 1,bd$myel_leukPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$myel_leukSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

haematological_cancer_function<-function(){

	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C81", "C82", "C83", "C84", "C85", "C86", "C87", "C88", "C89", "C90", "C91", "C92", "C93", "C94", "C95", "C96")
	haemICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(haemICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("200", "201", "202", "203", "204", "205", "206", "207", "208")
	haemICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	haemICD9 <-as.integer(haemICD9)
	table(haemICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = haemICD9,sitename = "haemICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = haemICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = haemICD10,sitename = "haemICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = haemICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$haemPrevelent<-ifelse(!is.na(bd$haemICD101), 1,NA)
	bd$haemPrevelent<-ifelse(!is.na(bd$haemICD91), 1,bd$haemPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$haemSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

head_and_neck_cancer<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C003", "C004", "C005", "C006", "C009", "C01", "C020", "C021", "C022", "C023", "C024", "C028", "C029", "C030", "C031", "C039", "C040", "C041", "C048", "C049", "C050", "C051", "C052", "C058", "C059", "C060", "C061", "C062", "C068", "C069", "C090", "C091", "C098", "C099", "C100", "C101", "C102", "C103", "C104", "C108", "C109", "C12", "C130", "C131", "C132", "C139", "C140", "C142", "C320", "C321", "C322", "C323", "C328", "C329")
	headneckICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(headneckICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("1410", "1412", "1413", "1419", "1430", "1431", "1449", "1450", "1451", "1452", "1453", "1455", "1460", "1461", "1610")
	headneckICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	headneckICD9 <-as.integer(headneckICD9)
	table(headneckICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = headneckICD9,sitename = "headneckICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = headneckICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = headneckICD10,sitename = "headneckICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = headneckICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$headneckPrevelent<-ifelse(!is.na(bd$headneckICD101), 1,NA)
	bd$headneckPrevelent<-ifelse(!is.na(bd$headneckICD91), 1,bd$headneckPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$headneckSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}



larynx_cancer<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C320", "C321", "C322", "C323", "C328", "C329")
	larynxICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(larynxICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("1610")
	larynxICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	larynxICD9 <-as.integer(larynxICD9)
	table(larynxICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = larynxICD9,sitename = "larynxICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = larynxICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = larynxICD10,sitename = "larynxICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = larynxICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$larynxPrevelent<-ifelse(!is.na(bd$larynxICD101), 1,NA)
	bd$larynxPrevelent<-ifelse(!is.na(bd$larynxICD91), 1,bd$larynxPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$larynxSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}


leukemia_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C91", "C92", "C93", "C94", "C95")
	leukICD10 <- unique (grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(leukICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("204", "205", "206", "207", "208")
	leukICD9 <- unique (grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	leukICD9 <-as.integer(leukICD9)
	table(leukICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = leukICD9,sitename = "leukICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = leukICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = leukICD10,sitename = "leukICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = leukICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$leukPrevelent<-ifelse(!is.na(bd$leukICD101), 1,NA)
	bd$leukPrevelent<-ifelse(!is.na(bd$leukICD91), 1,bd$leukPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$leukSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

liver_bile_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C22")
	liver_bileICD10 <- unique (grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(liver_bileICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("155")
	liver_bileICD9 <- unique (grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	liver_bileICD9 <-as.integer(liver_bileICD9)
	table(liver_bileICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = liver_bileICD9,sitename = "liver_bileICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = liver_bileICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = liver_bileICD10,sitename = "liver_bileICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = liver_bileICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$liver_bilePrevelent<-ifelse(!is.na(bd$liver_bileICD101), 1,NA)
	bd$liver_bilePrevelent<-ifelse(!is.na(bd$liver_bileICD91), 1,bd$liver_bilePrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$liver_bileSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

liver_cell_cancer_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C220")
	liver_cellICD10 <- unique (grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(liver_cellICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("1550")
	liver_cellICD9 <- unique (grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	liver_cellICD9 <-as.integer(liver_cellICD9)
	table(liver_cellICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = liver_cellICD9,sitename = "liver_cellICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = liver_cellICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = liver_cellICD10,sitename = "liver_cellICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = liver_cellICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$liver_cellPrevelent<-ifelse(!is.na(bd$liver_cellICD101), 1,NA)
	bd$liver_cellPrevelent<-ifelse(!is.na(bd$liver_cellICD91), 1,bd$liver_cellPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$liver_cellSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

lymphoid_leukemia_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C91")
	lymph_leukICD10 <- unique (grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(lymph_leukICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("204")
	lymph_leukICD9 <- unique (grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	lymph_leukICD9 <-as.integer(lymph_leukICD9)
	table(lymph_leukICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = lymph_leukICD9,sitename = "lymph_leukICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = lymph_leukICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = lymph_leukICD10,sitename = "lymph_leukICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = lymph_leukICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$lymph_leukPrevelent<-ifelse(!is.na(bd$lymph_leukICD101), 1,NA)
	bd$lymph_leukPrevelent<-ifelse(!is.na(bd$lymph_leukICD91), 1,bd$lymph_leukPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$lymph_leukSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

multiple_myeloma_function<-function(){
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C900")
	mult_myelICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(mult_myelICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c("2030")
	mult_myelICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	mult_myelICD9 <-as.integer(mult_myelICD9)
	table(mult_myelICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = mult_myelICD9,sitename = "mult_myelICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = mult_myelICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = mult_myelICD10,sitename = "mult_myelICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = mult_myelICD10,sitename = "otherICD10", other=T)

	print("functions complete!")

	# define overall
	bd$mult_myelPrevelent<-ifelse(!is.na(bd$mult_myelICD101), 1,NA)
	bd$mult_myelPrevelent<-ifelse(!is.na(bd$mult_myelICD91), 1,bd$mult_myelPrevelent)

	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)

	print("getting self-report cancers...")

	#self report
	#identify self reported all cancer (coded as 1)
	bd$mult_myelSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
												ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																		ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}

colorectal_cancer_function<-function(){
	#colorectal cancer
	allICD10 <- unlist(bd %>% select(starts_with("f.40006.")), use.names=F)
	allICD10<-allICD10[!is.na(allICD10)]
	allICD10<-allICD10[!duplicated(allICD10)]

	# subet for only C codes (D and O codes will be in "other" for exclusion from controls)
	allICD10 <- grep("C", allICD10, value=T)
	toMatch <- c("C180", "C181", "C182", "C183", "C184", "C185", "C186", "C187", "C188", "C189", "C19", "C20")
	colorectalICD10 <- unique(grep(paste(toMatch,collapse="|"), allICD10, value=TRUE))
	table(colorectalICD10)

	allICD9<-unlist(bd %>% select(starts_with("f.40013.")), use.names=F)
	allICD9<-allICD9[!is.na(allICD9)]
	allICD9<-allICD9[!duplicated(allICD9)]

	toMatch <- c(1530, 1531, 1532, 1533, 1534, 1535, 1536, 1537, 1538, 1539)
	colorectalICD9 <- unique(grep(paste(toMatch,collapse="|"), allICD9, value=TRUE))
	colorectalICD9 <-as.integer(colorectalICD9)
	table(colorectalICD9)

	bd<-UKBcancerFunc(dat=bd,cancerCode = colorectalICD9,sitename = "colorectalICD9", other=F, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = colorectalICD9,sitename = "otherICD9", other=T, cancer_col = "f.40013.")
	bd<-UKBcancerFunc(dat=bd,cancerCode = colorectalICD10,sitename = "colorectalICD10", other=F)
	bd<-UKBcancerFunc(dat=bd,cancerCode = colorectalICD10,sitename = "otherICD10", other=T)
	# define overall
	bd$colorectalPrevelent<-ifelse(!is.na(bd$colorectalICD101), 1,NA)
	bd$colorectalPrevelent<-ifelse(!is.na(bd$colorectalICD91), 1,bd$colorectalPrevelent)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD101), 1,NA)
	bd$otherPrevelent<-ifelse(!is.na(bd$otherICD91), 1,bd$otherPrevelent)
	bd$colorectalSelfreport<-ifelse(between(bd$f.20001.0.0, 1001, 99999),1,
							ifelse(between(bd$f.20001.0.1, 1001, 99999),1,
								ifelse(between(bd$f.20001.0.2, 1001, 99999),1,
									ifelse(between(bd$f.20001.0.3, 1001, 99999),1,
										ifelse(between(bd$f.20001.0.4, 1001, 99999),1,
											ifelse(between(bd$f.20001.0.5, 1001, 99999),1,
												ifelse(between(bd$f.20001.1.0, 1001, 99999),1,
													ifelse(between(bd$f.20001.1.1, 1001, 99999),1,
														ifelse(between(bd$f.20001.1.2, 1001, 99999),1,
															ifelse(between(bd$f.20001.1.3, 1001, 99999),1,
																ifelse(between(bd$f.20001.1.4, 1001, 99999),1,
																	ifelse(between(bd$f.20001.1.5, 1001, 99999),1,0))))))))))))
	return(bd)
}



format_smoking<-function(){
	lvl.0090 <- c(-3,0,1,2)
	lbl.0090 <- c("Prefer not to answer","Never","Previous","Current")
	bd$f.20116.0.0 <- ordered(bd$f.20116.0.0, levels=lvl.0090, labels=lbl.0090)
	bd$f.20116.1.0 <- ordered(bd$f.20116.1.0, levels=lvl.0090, labels=lbl.0090)
	bd$f.20116.2.0 <- ordered(bd$f.20116.2.0, levels=lvl.0090, labels=lbl.0090)
	names(bd)[names(bd) == "f.20116.0.0"]<-"smoking"
	return(bd)
}

format_sex<-function(){
	lvl.0009 <- c(0,1)
	lbl.0009 <- c("Female","Male")
	bd$f.31.0.0 <- as.character(ordered(bd$f.31.0.0, levels=lvl.0009, labels=lbl.0009))
	# names(bd)[names(bd)=="f.31.0.0"]<-"sex"
	return(bd)
	# names(bd)[names(bd) == "sex" ]<-"f.31.0.0"
}

format_behaviour<-function(){
	lvl.0039 <- c(-1,0,1,2,3,5,6,9)
	lbl.0039 <- c("Malignant","Benign","Uncertain whether benign or malignant","Carcinoma in situ","Malignant, primary site","Malignant, microinvasive","Malignant, metastatic site","Malignant, uncertain whether primary or metastatic site")	
	# bd$f.40012.0.0 <- ordered(bd$f.40012.0.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.1.0 <- ordered(bd$f.40012.1.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.2.0 <- ordered(bd$f.40012.2.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.3.0 <- ordered(bd$f.40012.3.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.4.0 <- ordered(bd$f.40012.4.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.5.0 <- ordered(bd$f.40012.5.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.6.0 <- ordered(bd$f.40012.6.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.7.0 <- ordered(bd$f.40012.7.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.8.0 <- ordered(bd$f.40012.8.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.9.0 <- ordered(bd$f.40012.9.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.10.0 <- ordered(bd$f.40012.10.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.11.0 <- ordered(bd$f.40012.11.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.12.0 <- ordered(bd$f.40012.12.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.13.0 <- ordered(bd$f.40012.13.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.14.0 <- ordered(bd$f.40012.14.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.15.0 <- ordered(bd$f.40012.15.0, levels=lvl.0039, labels=lbl.0039)
	# bd$f.40012.16.0 <- ordered(bd$f.40012.16.0, levels=lvl.0039, labels=lbl.0039)
	bd$f.40012.0.0 <- as.character(ordered(bd$f.40012.0.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.1.0 <- as.character(ordered(bd$f.40012.1.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.2.0 <- as.character(ordered(bd$f.40012.2.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.3.0 <- as.character(ordered(bd$f.40012.3.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.4.0 <- as.character(ordered(bd$f.40012.4.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.5.0 <- as.character(ordered(bd$f.40012.5.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.6.0 <- as.character(ordered(bd$f.40012.6.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.7.0 <- as.character(ordered(bd$f.40012.7.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.8.0 <- as.character(ordered(bd$f.40012.8.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.9.0 <- as.character(ordered(bd$f.40012.9.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.10.0 <- as.character(ordered(bd$f.40012.10.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.11.0 <- as.character(ordered(bd$f.40012.11.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.12.0 <- as.character(ordered(bd$f.40012.12.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.13.0 <- as.character(ordered(bd$f.40012.13.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.14.0 <- as.character(ordered(bd$f.40012.14.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.15.0 <- as.character(ordered(bd$f.40012.15.0, levels=lvl.0039, labels=lbl.0039))
	bd$f.40012.16.0 <- as.character(ordered(bd$f.40012.16.0, levels=lvl.0039, labels=lbl.0039))
	return(bd)
}

format_date_enrollment<-function(){
	bd$f.200.0.0 <- as.Date(bd$f.200.0.0)
	# names(bd)[names(bd)=="f.200.0.0"]<-"consenting_to_uk"
	return(bd)	
}

format_date_of_attending_assessment_centre<-function(){
	bd$f.53.0.0 <- as.Date(bd$f.53.0.0)
	bd$f.53.1.0 <- as.Date(bd$f.53.1.0)
	bd$f.53.2.0 <- as.Date(bd$f.53.2.0)
	bd$f.53.3.0 <- as.Date(bd$f.53.3.0)
	return(bd)
}

format_date_of_death<-function(){
	bd$f.40000.0.0 <- as.Date(bd$f.40000.0.0)
	bd$f.40000.1.0 <- as.Date(bd$f.40000.1.0)
	names(bd)[names(bd)=="f.40000.0.0"]<-"date_of_death_40000"
	return(bd)
	
}

# unique(bd$f.40005.0.0)
format_date_diagnosis<-function(){
	# bd$f.40005.0.0
	bd$f.40005.0.0 <- as.Date(bd$f.40005.0.0)
	bd$f.40005.1.0 <- as.Date(bd$f.40005.1.0)
	bd$f.40005.2.0 <- as.Date(bd$f.40005.2.0)
	bd$f.40005.3.0 <- as.Date(bd$f.40005.3.0)
	bd$f.40005.4.0 <- as.Date(bd$f.40005.4.0)
	bd$f.40005.5.0 <- as.Date(bd$f.40005.5.0)
	bd$f.40005.6.0 <- as.Date(bd$f.40005.6.0)
	bd$f.40005.7.0 <- as.Date(bd$f.40005.7.0)
	bd$f.40005.8.0 <- as.Date(bd$f.40005.8.0)
	bd$f.40005.9.0 <- as.Date(bd$f.40005.9.0)
	bd$f.40005.10.0 <- as.Date(bd$f.40005.10.0)
	bd$f.40005.11.0 <- as.Date(bd$f.40005.11.0)
	bd$f.40005.12.0 <- as.Date(bd$f.40005.12.0)
	bd$f.40005.13.0 <- as.Date(bd$f.40005.13.0)
	bd$f.40005.14.0 <- as.Date(bd$f.40005.14.0)
	bd$f.40005.15.0 <- as.Date(bd$f.40005.15.0)
	bd$f.40005.16.0 <- as.Date(bd$f.40005.16.0)
	return(bd)
}


genotyping_batch<-function(){
	names(bd)[names(bd)=="f.22000.0.0"]<-"genotyping_batch"
	return(bd)
}

format_nsaids_baseline_6154<-function(){
	# Data-Field 6154
	# Medication for pain relief, constipation, heartburn
	# Do you regularly take any of the following? (You can select more than one answer)
	lvl.100628 <- c(-7,-3,-1,1,2,3,4,5,6)
	lbl.100628 <- c("None of the above","Prefer not to answer","Do not know","Aspirin","Ibuprofen (e.g. Nurofen)","Paracetamol","Ranitidine (e.g. Zantac)","Omeprazole (e.g. Zanprol)","Laxatives (e.g. Dulcolax, Senokot)")
	bd$f.6154.0.0 <- ordered(bd$f.6154.0.0, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.0.1 <- ordered(bd$f.6154.0.1, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.0.2 <- ordered(bd$f.6154.0.2, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.0.3 <- ordered(bd$f.6154.0.3, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.0.4 <- ordered(bd$f.6154.0.4, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.0.5 <- ordered(bd$f.6154.0.5, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.1.0 <- ordered(bd$f.6154.1.0, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.1.1 <- ordered(bd$f.6154.1.1, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.1.2 <- ordered(bd$f.6154.1.2, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.1.3 <- ordered(bd$f.6154.1.3, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.1.4 <- ordered(bd$f.6154.1.4, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.1.5 <- ordered(bd$f.6154.1.5, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.2.0 <- ordered(bd$f.6154.2.0, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.2.1 <- ordered(bd$f.6154.2.1, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.2.2 <- ordered(bd$f.6154.2.2, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.2.3 <- ordered(bd$f.6154.2.3, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.2.4 <- ordered(bd$f.6154.2.4, levels=lvl.100628, labels=lbl.100628)
	bd$f.6154.2.5 <- ordered(bd$f.6154.2.5, levels=lvl.100628, labels=lbl.100628)

	# any NSAIDS
	Test_list<-NULL
	for(i in 1:nrow(bd)){
	# for(i in 1:100){
		print(i)
		Temp<-bd[i,c("f.6154.0.0","f.6154.0.1","f.6154.0.2","f.6154.0.3","f.6154.0.4","f.6154.0.5")]
		Test_list[[i]]<-any(Temp == "Aspirin" | Temp == "Ibuprofen (e.g. Nurofen)")
	}

	Test_list<-unlist(Test_list)
	Test_list[is.na(Test_list)]<-"no"
	Test_list[Test_list=="TRUE"]<-"yes"
	bd$nsaid_baseline_f_6154<-NA
	bd$nsaid_baseline_f_6154[Test_list=="no"]<-"no"
	bd$nsaid_baseline_f_6154[Test_list=="yes"]<-"yes"
	# unique(bd$nsaid_baseline_f_6154)
	return(bd)
}
	# table(bd1$nsaid_baseline_f_6154,bd1$nsaid_baseline_f_20003)
	# table(bd1$aspirin_baseline_f_6154,bd1$aspirin_baseline_f_20003)

format_aspirin_baseline_6154<-function(){
	#Aspirin
	Test_list<-NULL
	for(i in 1:nrow(bd)){
		print(i)
		Temp<-bd[i,c("f.6154.0.0","f.6154.0.1","f.6154.0.2","f.6154.0.3","f.6154.0.4","f.6154.0.5")]
		Test_list[[i]]<-any(Temp == "Aspirin")
	}

	Test_list2<-unlist(Test_list)
	Test_list2[is.na(Test_list2)]<-"no"
	Test_list2[Test_list2=="TRUE"]<-"yes"
	bd$aspirin_baseline_f_6154<-NA
	bd$aspirin_baseline_f_6154[Test_list2=="no"]<-"no"
	bd$aspirin_baseline_f_6154[Test_list2=="yes"]<-"yes"

	# pilot study
	lvl.100688 <- c(-7,-3,-1,1,2,3,4,5)
	lbl.100688 <- c("None of the above","Prefer not to answer","Do not know","Aspirin","Ibuprofen (e.g. Nurofen)","Paracetamol","Codeine","Ranitidine (e.g. Zantac)")
	bd$f.10004.0.0 <- ordered(bd$f.10004.0.0, levels=lvl.100688, labels=lbl.100688)
	bd$f.10004.0.1 <- ordered(bd$f.10004.0.1, levels=lvl.100688, labels=lbl.100688)
	bd$f.10004.0.2 <- ordered(bd$f.10004.0.2, levels=lvl.100688, labels=lbl.100688)
	bd$f.10004.0.3 <- ordered(bd$f.10004.0.3, levels=lvl.100688, labels=lbl.100688)
	bd$f.10004.0.4 <- ordered(bd$f.10004.0.4, levels=lvl.100688, labels=lbl.100688)
	return(bd)
}

med_codings<-function(){
	Med<-readLines("~/UKBB_cancer_outcomes/coding4.tsv")
	Med2<-data.frame(do.call(rbind,strsplit(Med,split="\t")))
	names(Med2)<-paste(Med2[1,])
	Med2<-Med2[2:nrow(Med2),]
	nsaids<-c("ibuprofen","naproxen","diclofenac","celecoxib","mefenamic acid","etoricoxib","indomethacin","aspirin") #https://www.nhs.uk/conditions/nsaids/
	Pos<-unique(unlist(lapply(nsaids,FUN=function(x) grep(x,Med2$meaning))))
	Med2<-Med2[Pos,]
	return(Med2)
}

format_nsaids_baseline_20003<-function(){

# f.20003
# Treatment/medication code
# verbal interview
	# This category contains data obtained through a verbal interview by a trained nurse on prescription medications and includes data on type and number of medications taken.
	# The interviewer was made aware, via a pop-up box on their computer screen, if the participant had answered in the touchscreen that they are taking regular prescription medication, and was then prompted to ask "Could you now tell me what these are?" If the participant indicated in the touchscreen that they were taking any of the following classes of medications: blood pressure lowering, cholesterol lowering, hormone replacement therapy or oral contraceptive pills, then the interviewer was prompted to record the name of the medication. If the participant stated in the touchscreen they were not taking any regular prescription medications (or were not sure), this question was asked again and confirmed by the interviewer.
	# This category contains data on any regular treatments taken weekly, monthly, etc. It does not include short-term medications (such as a 1 week course of antibiotics) or prescribed medication that is not taken, or over-the-counter medications, vitamins and supplements (this information was collected in the touchscreen and was not recorded here, unless for some reason the participant had forgotten to record it in the touchscreen). Doses and formulations were not recorded.
	# Medicines that could not be coded at the time of the interview were entered as free text, and subsequently coded wherever possible.
	
	# old code
	# nsaid_list<-NULL
	# length(which(!is.na(bd1$f.20003.0.3)))
	# for(i in 1:length(Names)){
	# 	print(i)
	# 	print(Names[i])
	# 	bd1<-merge(bd,Med2,by.x=Names[i],by.y="coding",all.x=T)
	# 	Col<-paste0("nsaid",i)
	# 	bd1[,Col]<-NA		
	# 	bd1[,Col][is.na(bd1$meaning)]<-"no"
	# 	bd1[,Col][!is.na(bd1$meaning)]<-"yes"
	# 	nsaid_list[[i]]<-bd1[,Col]
	# }
	# nsaids_list2<-data.frame(do.call(cbind,nsaid_list))
	# Test_any_yes<-NULL
	# # Test_all_yes<-NULL
	# for(i in 1:nrow(nsaids_list2)){
	# # for(i in 1:1000){
	# 	print(i)
	# 	Test_any_yes[[i]]<-any(nsaids_list2[i,]=="yes")
	# 	# Test_all_yes[[i]]<-all(nsaids_list2[i,]=="yes")
	# }

	# Test_any_yes2<-unlist(Test_any_yes)	
	# bd1$nsaid_baseline_f_20003<-NA
	# bd1$nsaid_baseline_f_20003[!Test_any_yes2]<-"no"
	# bd1$nsaid_baseline_f_20003[Test_any_yes2]<-"yes"

	Med2<-med_codings()
	nsaid_codings<-as.numeric(Med2$coding)

	Names<-names(bd)[grep("20003.0",names(bd))] 
	#####################
	# Any nsaids f.20003#
	#####################
	Test_list<-NULL
	for(i in 1:nrow(bd)){
	# for(i in 1:100){
		# i<-1
		print(i) 
		Temp<-bd[,Names]
		Test_list[[i]]<-any(Temp[i,] %in% c(nsaid_codings))
	}

	Test_list2<-unlist(Test_list)
	bd$nsaid_baseline_f_20003<-NA
	bd$nsaid_baseline_f_20003[!Test_list2]<-"no"
	bd$nsaid_baseline_f_20003[Test_list2]<-"yes"
	return(bd)
}

format_aspirin_baseline_20003<-function(){
	Med2<-med_codings()
	aspirin_codings<-as.numeric(Med2$coding[grep("aspirin",Med2$meaning)])
	Names<-names(bd)[grep("20003.0",names(bd))] 
	#####################
	# aspirin f.20003#
	#####################
	Test_list<-NULL
	for(i in 1:nrow(bd)){
	# for(i in 1:100){
		print(i)
		Temp<-bd[,Names]
		Test_list[[i]]<-any(Temp[i,] %in% c(aspirin_codings))
	}
	Test_list2<-unlist(Test_list)
	bd$aspirin_baseline_f_20003<-NA
	bd$aspirin_baseline_f_20003[!Test_list2]<-"no"
	bd$aspirin_baseline_f_20003[Test_list2]<-"yes"	
	return(bd)
}

cleanup_names<-function(){
	
	bd <- df_split %>% select("projectID", "geneticID", starts_with("f."), starts_with("incident"),starts_with("overall"))
	names(bd)
	bd<-bd[,!names(bd) %in% c("incident.flag","overall_cases", "overall_cases2"  )]

	# Names_keep<-c(names(df_split)[grep("f\\.",names(df_split))],"projectID","geneticID","allSelfreport"  ,"incident_pan_inclc44_cancer","overall_pan_inclc44_cancer" )

	# Names_keep<-c(names(df_split)[grep("f\\.",names(df_split))],"projectID","geneticID","allSelfreport"  ,starts_with("incident"),starts_with("overall"))
	# bd<-df_split[,names(df_split) %in% Names_keep]
	return(bd)
}

lung_cancer_function2<-function(dat=NULL){
	setwd("~/UKBB_cancer_outcomes")
	######################################################### Format results #################################################################################

	# 1. split the cancer diagnoses: #ICD_code/date_of_diagnosis/histology_code/behaviour_code/age_at_diagnosis
	# 2. format columns
	# 3. generate incidence of cancer flag
	# 4. generate tumour behaviour flag 
	# 5. define controls
	# 6. define incident cases
	# 7. define overall cases
	# 8. tidy up data

	library(tidyr); library(dplyr)
	#1. separate the cancer data into columns

	print("performing task 1. separation...")
	Df<-dat
	Df2 <-separate(Df, lungICD91, into = c("lung.ICD9.1", "lung.ICD9.date.diagnosis.1", "lung.ICD9.histology.1", "lung.ICD9.behaviour.1", "lung.ICD9.age_diagnosis.1"), sep = "/")
	Df3 <-separate(Df2, lungICD92, into = c("lung.ICD9.2", "lung.ICD9.date.diagnosis.2", "lung.ICD9.histology.2", "lung.ICD9.behaviour.2", "lung.ICD9.age_diagnosis.2"), sep = "/")

	Df4 <-separate(Df3, lungICD101, into = c("lung.ICD10.1", "lung.ICD10.date.diagnosis.1", "lung.ICD10.histology.1", "lung.ICD10.behaviour.1", "lung.ICD10.age_diagnosis.1"), sep = "/")
	Df5 <-separate(Df4, lungICD102, into = c("lung.ICD10.2", "lung.ICD10.date.diagnosis.2", "lung.ICD10.histology.2", "lung.ICD10.behaviour.2", "lung.ICD10.age_diagnosis.2"), sep = "/")
	Df6 <-separate(Df5, lungICD103, into = c("lung.ICD10.3", "lung.ICD10.date.diagnosis.3", "lung.ICD10.histology.3", "lung.ICD10.behaviour.3", "lung.ICD10.age_diagnosis.3"), sep = "/")

	Df7 <-separate(Df6, otherICD91, into = c("other.ICD9.1", "other.ICD9.date.diagnosis.1", "other.ICD9.histology.1", "other.ICD9.behaviour.1", "other.ICD9.age_diagnosis.1"), sep = "/")
	Df8 <-separate(Df7, otherICD92, into = c("other.ICD9.2", "other.ICD9.date.diagnosis.2", "other.ICD9.histology.2", "other.ICD9.behaviour.2", "other.ICD9.age_diagnosis.2"), sep = "/")
	Df9 <-separate(Df8, otherICD93, into = c("other.ICD9.3", "other.ICD9.date.diagnosis.3", "other.ICD9.histology.3", "other.ICD9.behaviour.3", "other.ICD9.age_diagnosis.3"), sep = "/")
	Df10 <-separate(Df9, otherICD94, into = c("other.ICD9.4", "other.ICD9.date.diagnosis.4", "other.ICD9.histology.4", "other.ICD9.behaviour.4", "other.ICD9.age_diagnosis.4"), sep = "/")
	Df11 <-separate(Df10, otherICD95, into = c("other.ICD9.5", "other.ICD9.date.diagnosis.5", "other.ICD9.histology.5", "other.ICD9.behaviour.5", "other.ICD9.age_diagnosis.5"), sep = "/")
	Df12 <-separate(Df11, otherICD96, into = c("other.ICD9.6", "other.ICD9.date.diagnosis.6", "other.ICD9.histology.6", "other.ICD9.behaviour.6", "other.ICD9.age_diagnosis.6"), sep = "/")
	Df13 <-separate(Df12, otherICD97, into = c("other.ICD9.7", "other.ICD9.date.diagnosis.7", "other.ICD9.histology.7", "other.ICD9.behaviour.7", "other.ICD9.age_diagnosis.7"), sep = "/")
	Df14 <-separate(Df13, otherICD98, into = c("other.ICD9.8", "other.ICD9.date.diagnosis.8", "other.ICD9.histology.8", "other.ICD9.behaviour.8", "other.ICD9.age_diagnosis.8"), sep = "/")

	Df15 <-separate(Df14, otherICD101, into = c("other.ICD10.1", "other.ICD10.date.diagnosis.1", "other.ICD10.histology.1", "other.ICD10.behaviour.1", "other.ICD10.age_diagnosis.1"), sep = "/")
	Df16 <-separate(Df15, otherICD102, into = c("other.ICD10.2", "other.ICD10.date.diagnosis.2", "other.ICD10.histology.2", "other.ICD10.behaviour.2", "other.ICD10.age_diagnosis.2"), sep = "/")
	Df17 <-separate(Df16, otherICD103, into = c("other.ICD10.3", "other.ICD10.date.diagnosis.3", "other.ICD10.histology.3", "other.ICD10.behaviour.3", "other.ICD10.age_diagnosis.3"), sep = "/")
	Df18 <-separate(Df17, otherICD104, into = c("other.ICD10.4", "other.ICD10.date.diagnosis.4", "other.ICD10.histology.4", "other.ICD10.behaviour.4", "other.ICD10.age_diagnosis.4"), sep = "/")
	Df19 <-separate(Df18, otherICD105, into = c("other.ICD10.5", "other.ICD10.date.diagnosis.5", "other.ICD10.histology.5", "other.ICD10.behaviour.5", "other.ICD10.age_diagnosis.5"), sep = "/")
	Df20 <-separate(Df19, otherICD106, into = c("other.ICD10.6", "other.ICD10.date.diagnosis.6", "other.ICD10.histology.6", "other.ICD10.behaviour.6", "other.ICD10.age_diagnosis.6"), sep = "/")
	df_split <-separate(Df20, otherICD107, into = c("other.ICD10.7", "other.ICD10.date.diagnosis.7", "other.ICD10.histology.7", "other.ICD10.behaviour.7", "other.ICD10.age_diagnosis.7"), sep = "/")

	rm(list = c("Df2", "Df3", "Df4", "Df5", "Df6", "Df7", "Df8", "Df9", "Df10", "Df11", "Df12", "Df13", "Df14", "Df15", "Df16", "Df17", "Df18", "Df19", "Df20"))

	str(Df, list.len=ncol(Df))
	str(df_split, list.len=ncol(df_split))

	# rm(Df)

	# 2. format the columns

	print("2. formatting columns...")

	df_split$enroll <-as.Date(df_split$f.200.0.0)
	date <- grepl("date", names(df_split))
	df_split[,date] <-lapply(df_split[, date, drop=FALSE], as.Date)
	age <- grepl("age", names(df_split))
	df_split[,age] <-lapply(df_split[, age, drop=FALSE], as.numeric)

	# 3. generate incident cancer flags: incidence is classed as cancer cases diagnosed after enrolment to UKBB (var: f.200.0.0, Date of consenting to join UK Biobank)
	print("generating incidence flag")
	# Get earliest date for the cancer (only need the first instance)
	df_split <- df_split %>% mutate(earliest_date = pmin(lung.ICD9.date.diagnosis.1, lung.ICD10.date.diagnosis.1, na.rm =T))

	df_split$incident.flag <- ifelse(df_split$earliest_date >= df_split$enroll,1,0)   
	table(df_split$incident.flag)

	# 4. generate behaviour flag: only using codes: 3, 6, 7, & 9 see below
	# behaviour levels: "Malignant, primary site","Malignant, microinvasive","Malignant, metastatic site","Malignant, uncertain whether primary or metastatic site"
	# restrict to "Malignant, primary site"


	# vast majority of cancers are "Malignant, primary site"
	table(df_split$lung.ICD9.behaviour.1)
	table(df_split$lung.ICD9.behaviour.2)
	table(df_split$lung.ICD10.behaviour.1)
	table(df_split$lung.ICD10.behaviour.2)
	table(df_split$lung.ICD10.behaviour.3)
 
	print("generating behaviour flag")
	df_split$behaviour.flag <-
	  ifelse(
	    df_split$lung.ICD9.behaviour.1 == "Malignant, primary site" |
	      df_split$lung.ICD9.behaviour.2 == "Malignant, primary site" |
	      df_split$lung.ICD10.behaviour.1 == "Malignant, primary site" |
	      df_split$lung.ICD10.behaviour.2 == "Malignant, primary site" |
	      df_split$lung.ICD10.behaviour.3 == "Malignant, primary site"
	      ,1,0)

	# df_split$behaviour.flag <-
	#   ifelse(
	#     df_split$lung.ICD9.behaviour.1 == "Malignant, primary site" |
	#       df_split$lung.ICD9.behaviour.1 == "Malignant, microinvasive" |
	#       df_split$lung.ICD9.behaviour.1 == "Malignant, metastatic site" |
	#       df_split$lung.ICD9.behaviour.1 == "Malignant, uncertain whether primary or metastatic site" |
	#       df_split$lung.ICD9.behaviour.2 == "Malignant, primary site" |
	#       df_split$lung.ICD9.behaviour.2 == "Malignant, microinvasive" |
	#       df_split$lung.ICD9.behaviour.2 == "Malignant, metastatic site" |
	#       df_split$lung.ICD9.behaviour.2 == "Malignant, uncertain whether primary or metastatic site"    |
	#       df_split$lung.ICD10.behaviour.1 == "Malignant, primary site" |
	#       df_split$lung.ICD10.behaviour.1 == "Malignant, microinvasive" |
	#       df_split$lung.ICD10.behaviour.1 == "Malignant, metastatic site" |
	#       df_split$lung.ICD10.behaviour.1 == "Malignant, uncertain whether primary or metastatic site"     |
	#       df_split$lung.ICD10.behaviour.2 == "Malignant, primary site" |
	#       df_split$lung.ICD10.behaviour.2 == "Malignant, microinvasive" |
	#       df_split$lung.ICD10.behaviour.2 == "Malignant, metastatic site" |
	#       df_split$lung.ICD10.behaviour.2 == "Malignant, uncertain whether primary or metastatic site"   |
	#       df_split$lung.ICD10.behaviour.3 == "Malignant, primary site" |
	#       df_split$lung.ICD10.behaviour.3 == "Malignant, microinvasive" |
	#       df_split$lung.ICD10.behaviour.3 == "Malignant, metastatic site" |
	#       df_split$lung.ICD10.behaviour.3 == "Malignant, uncertain whether primary or metastatic site",
	#     1,
	#     0
	#   )

	table(df_split$behaviour.flag)

	# 5. generate controls: controls are participants that do not have a cancer of interest code or any other cancer code including ICD10:D codes
	# controls also have no self-report of cancers

	print("defining controls...")
	# control flags (controls =1, others =0)
	df_split$controls <- ifelse(
	  is.na(df_split$lung.ICD9.1) &
	    is.na(df_split$lung.ICD9.2) &
	    is.na(df_split$other.ICD9.1) &
	    is.na(df_split$other.ICD9.2) &
	    is.na(df_split$other.ICD9.3) &
	    is.na(df_split$other.ICD9.4) &
	    is.na(df_split$other.ICD9.5) &
	    is.na(df_split$other.ICD9.6) &
	    is.na(df_split$other.ICD9.7) &
	    is.na(df_split$other.ICD9.8) &
	    is.na(df_split$lung.ICD10.1) &
	    is.na(df_split$lung.ICD10.2) &
	    is.na(df_split$lung.ICD10.3) &
	    is.na(df_split$other.ICD10.1) &
	    is.na(df_split$other.ICD10.2) &
	    is.na(df_split$other.ICD10.3) &
	    is.na(df_split$other.ICD10.4) &
	    is.na(df_split$other.ICD10.5) &
	    is.na(df_split$other.ICD10.6) &
	    is.na(df_split$other.ICD10.7) &
	    is.na(df_split$lungSelfreport) &
	    is.na(df_split$lungPrevelent) &
	    is.na(df_split$otherPrevelent),
	  1,
	  0
	)

	names(df_split)

	table(df_split$controls)


	# 5. generate incident cases: participants who have a cancer of interest code diagnosed after enrolment
	# Self report not included (some report NA but have a diagnosis, and some report cancer but these are carcinoma in situ ICD10:D codes)

	print("defining incident cases...")

	# define incident cases (2=cases)
	df_split$cases <- ifelse(df_split$incident.flag ==1 & df_split$behaviour.flag ==1, 2, 0) # cases

	table(df_split$cases)

	# make incident lung cancer outcome
	df_split$incident_lung_cancer <- NA
	df_split$incident_lung_cancer <- ifelse(df_split$controls==1, 1, ifelse(df_split$cases==2, 2, NA))

	table(df_split$incident_lung_cancer)

	# 6. make overall cancer outcome (these are cases diagnosed both before and after enrolment)
	print("generating overall cancer cases")

	# define overall cases (2=cases)
	#df_split$overall_cases <- ifelse(df_split$lungPrevelent ==1 & df_split$behaviour.flag ==1, 2, 0) # cases
	#table(df_split$overall_cases); table(df_split$overall_cases2)

	df_split$overall_cases <- ifelse(
	  !is.na(df_split$lung.ICD9.1) |
	    !is.na(df_split$lung.ICD9.2) |
	    !is.na(df_split$lung.ICD10.1) |
	    !is.na(df_split$lung.ICD10.2) |
	    !is.na(df_split$lung.ICD10.3),
	  2,
	  0
	)

	df_split$overall_cases2 <- ifelse(df_split$overall_cases ==2 & df_split$behaviour.flag ==1, 2, 0) # cases

	table(df_split$overall_cases); table(df_split$overall_cases2)

	# make overall lung cancer outcome
	df_split$overall_lung_cancer <- NA
	df_split$overall_lung_cancer <- ifelse(df_split$controls==1, 1, ifelse(df_split$overall_cases2==2, 2, NA))
	return(df_split)

	# table(df_split$overall_lung_cancer)

	# # 7. Numbers and tidying
	# bc <- df_split[,c("projectID", "f.31.0.0", "incident_lung_cancer", "overall_lung_cancer")]

	# #cases =2, controls =1
	# print("numbers for incident cancer")
	# table(bc$incident_lung_cancer)

	# print("numbers for overall cancer")
	# table(bc$overall_lung_cancer)

	# # formatting for BOLT LMM pipeline:
	# # a. link with genetic IEU IDs
	# library(readr)
	# linker <- read_csv("../linker.csv") #cols = ieu, app
	# print("linking IDs")
	# bc <- merge(bc, linker, by.x = "projectID", by.y = "app") #not all match re-do the numbers

	# #cases =2, controls =1
	# print("numbers for incident cancer")
	# table(bc$incident_lung_cancer)

	# print("numbers for overall cancer")
	# table(bc$overall_lung_cancer)

	# #These files should be space delimited text files.
	# #• The first two columns must be FID and IID (the PLINK identifiers of an individual); any number of columns may follow.
	# #• Values in the column should be numeric.
	# #• Case/control phenotypes should be encoded as 1=unaffected (control), 2=affected (case).

	# bc$FID <-bc$ieu
	# bc$IID <-bc$ieu

	# ## Now merge with genetic samples and covariates (complete cases) to get actual case/control numbers for the GWAS
	# sample <- read.table("../sample.txt", header=T, stringsAsFactors=F)
	# covars <- read.table("../covariates.txt", header=T, stringsAsFactors=F)

	# df <- merge(bc, sample, by.x = "FID", by.y = "FID")
	# df <- merge(df, covars, by.x = "FID", by.y = "FID")
	# inc_df <- subset(df, df$incident_lung_cancer !="NA" & sex.y !="NA")
	# overall_df <- subset(df, df$overall_lung_cancer !="NA" & sex.y !="NA")


	# print("numbers for incident cancer")
	# table(inc_df$incident_lung_cancer)

	# print("numbers for overall cancer")
	# table(overall_df$overall_lung_cancer)

	# write.table(bc[,c("FID", "IID", "incident_lung_cancer", "overall_lung_cancer")], file="../UKBB_lung_cancer.txt", sep=" ", row.names = F, quote = F)
}