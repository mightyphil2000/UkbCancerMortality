#' Multiple cancer diagnoses 
#'
#' retrieve individuals with multiple occurences in the cancer registry (presumably corresponding to multiple cancer diagnoses), including for each cancer occurence: the date of diagnosis, the icd10 or icd9 code, the date of death, primary causes of death and secondary causes of death. The function retrieves a single example for each occurence of multiple diagnoses. The idea is to create examples where an individual occurs in the cancer registry multiple times so as to understand the data better. 
#'
#' @param dat the cancer dataset 
#' @param cancer_name the name of the field/column in dat corresponding to case-control status for the cancer of interest. The field must take on values of 1 (control) or 2 (cases). 
#' @param occurences the number of cancer diagnoses. Can be a single value or range of values, e.g. 2:5 will extract individuals with 2 to 5 cancer diagnoses. If set to 5, the function will retrieve one individual with at least 5 occurences in the cancer registry. If 2:5. then will retrieve one individual with 2, 3, 4 and 5 occurences in the cancer registry.   
#'
#' @return data frame
#' @export

multiple_diagnoses<-function(dat=NULL,cancer_name=NULL,occurences=NULL){
	dat<-dat[which(dat[,cancer_name] == 2),]	
	Names<-names(dat)[grep("f.40005",names(dat))]
	Names2<-names(dat)[grep("f.40006",names(dat))]
	Names3<-names(dat)[grep("f.40013",names(dat))]
	Names4<-names(dat)[grep("f.40001",names(dat))]
	Names5<-names(dat)[grep("f.40002",names(dat))]
	Names6<-names(dat)[grep("f.40000",names(dat))]
	N_diag<-lapply(1:nrow(dat),FUN=function(x) 
			unique(as.numeric(dat[x,Names])))
	dat$N_diagnoses<-unlist(lapply(1:nrow(dat),FUN=function(x) 
			length(unlist(N_diag[x])[!is.na(unlist(N_diag[x]))])))
	Test_list<-NULL
	for(i in occurences){
		Test<-dat[which(dat$N_diagnoses==i)[1],c("geneticID",Names,Names2,Names4,Names5,Names6)]
		Pos<-unlist(lapply(1:ncol(Test),FUN=function(i) !is.na(Test[,i])))
		Test[,Pos]
		Test_list[[i]]<-Test[,Pos]
	}
	dat2<-do.call(plyr::rbind.fill,Test_list)
	return(dat2)
}

# dat=lun3
multiple_diagnoses2<-function(dat=NULL,cancer_name=NULL,icd9=NULL,icd10=NULL){
	if(!is.null(cancer_name))
	{
		dat<-dat[which(dat[,cancer_name] == 2),] #a redundant step if dataset already restricted to target cancer	
	}
	Names<-names(dat)[grep("f.40005",names(dat))] #Date of diagnosis
	Names2<-names(dat)[grep("f.40006",names(dat))] #Cancer diagnosis ICD10
	Names3<-names(dat)[grep("f.40013",names(dat))] #Cancer diagnosis ICD9
	# Names4<-names(dat)[grep("f.40001",names(dat))] #Primary cause of death
	# Names5<-names(dat)[grep("f.40002",names(dat))] #Secondary cause of death
	# Names6<-names(dat)[grep("f.40000",names(dat))] #Date of death
	N_diag<-lapply(1:nrow(dat),FUN=function(x) 
			unique(as.numeric(dat[x,Names])))
	dat$N_diagnoses<-unlist(lapply(1:nrow(dat),FUN=function(x) 
			length(unlist(N_diag[x])[!is.na(unlist(N_diag[x]))])))
	dat1<-dat[which(dat$N_diagnoses>1),]
	dat1diag<-dat[which(dat$N_diagnoses==1),]
	dat1$date_diagnosis2 # this is the earliest diagnosis date for target cancer
	#first_diagnosis: this is the first/earliest diagnosis of any cancer. may or may not be the same as date_diagnosis2 if the first diagnosed cancer was different to the target cancer

	min_date_diagnosis<-unlist(lapply(1:nrow(dat1),FUN=function(x) 
			min(as.numeric(dat1[x,Names]),na.rm=TRUE)))
	dat1$first_diagnosis<-as.Date(min_date_diagnosis,origin="1970-01-01")
	# calculate difference in time between any cancer first diagnosis and target cancer first diagnosis. time is in days
	
	# where difference = 0 we assume target cancer was first cancer diagnosed (but need to check for other cancers diagnosed on same day)
	Index<-which(dat1$date_diagnosis2-dat1$first_diagnosis==0) 
	dat1$target_cancer_first<-0
	dat1$target_cancer_first[Index]<-1
	
	# identify second earliest cancer diagnosis date and check to see if this was <3 months of first diagnosis. 
	# which instance of f.40005 contains the earliest diagnosis. It is not always the first instance
	# Pos<-which(Instance_min_date==2)
	Instance_min_date<-lapply(1:nrow(dat1),FUN=function(x) 
		which(as.numeric(dat1[x,Names]) == min_date_diagnosis[x]))
	Pos<-unlist(lapply(1:length(Instance_min_date),FUN=function(x)
		length(unlist(Instance_min_date[x]))>1))
	
	dat2<-dat1[Pos,] #dataset where there is more than one diagnosis on the date of the earliest diagnosis (ie more than one diagnosis on same day)
	dat3<-dat1[!Pos,] #dataset where there is only one diagnosis on the date of the earliest diagnosis
	min_date_diagnosis_dat3<-min_date_diagnosis[!Pos]
	
	# min_date_diagnosis2_dat3<-NULL
	# Instance_min_dat<-NULL
	# Instance_min_date2<-NULL
	# ICD_list<-NULL
	# date_diag_list<-NULL
	# diff_date_list<-NULL
	

	Index_6mo<-NULL
	Index_3yr<-NULL
	for(i in 1:nrow(dat1)){
		# print(i)
		# Instance_min_date[[i]]<-which(as.numeric(dat1[i,Names]) == min_date_diagnosis[i])
		date_diag<-as.numeric(dat1[i,Names])
		# date_diag<-date_diag[!is.na(date_diag)]
		diff_date<-date_diag-as.numeric(dat1$date_diagnosis2[i])
		diff_date<-diff_date[!is.na(diff_date)]
		Pos<-which(!is.na(as.numeric(dat1[i,Names])))
		ICD10<-dat1[i,Names2]
		Pos10<-which(!is.na(ICD10))
		# ICD10<-ICD10[Pos10]
		ICD9<-dat1[i,Names3]
		Pos9<-which(!is.na(ICD9))
		# if(any(!is.na(ICD10)) & any(!is.na(ICD9)))
		# {
		# 	if(Pos9>Pos10) stop("instance of ICD9 is not before instance of icd10")
		# }
		ICD<-c(ICD9[Pos9],ICD10[Pos10])
		ICD<-unlist(ICD)
		ORDER<-c(Pos9,Pos10) 
		ICD<-ICD[order(ORDER)]
		days_six_months<-365.25/12*6
		days_three_years<-365.25*3

		Index_6mo[[i]]<-FALSE		
		# find non-target-cancer diagnoses that occured during the 6 month period after the target-cancer diagnosis
		if(any(unlist(diff_date) >=0 & unlist(diff_date) <  days_six_months ) )
		{
			Pos<-which(unlist(diff_date) >=0 & unlist(diff_date) <  days_six_months )
			# if(!all(unlist(diff_date) >=0 & unlist(diff_date) <  days_six_months)) stop("")
			
			if(length(Pos)>5) stop("more than three diagnosis during 6 month period after target cancer diagnosis")
			
			#Check that the second cancer diagnosis is not the non-target cancer
			if(any(icd9!=substring(ICD[Pos],1,2) & icd10!=substring(ICD[Pos],1,3)))
			{
				Index_6mo[[i]]<-TRUE
			}
			
		}

		Index_3yr[[i]]<-FALSE
		if(any(unlist(diff_date) <0 & unlist(diff_date) >  -1*days_three_years))
		{


			Pos<-which(unlist(diff_date) <0 & unlist(diff_date) >  -1*days_three_years)	
			if(length(Pos)>3) stop("more than three diagnoses in 3 period before target cancer diagnosis")

	
			if(any(icd9!=substring(ICD[Pos],1,2) & icd10!=substring(ICD[Pos],1,3)))
			{
				Index_3yr[[i]]<-TRUE	
			}

		}
	}

	dat1$Index_6mo<-unlist(Index_6mo)
	dat1$Index_3yr<-unlist(Index_3yr)
	dat1diag$Index_6mo<-FALSE
	dat1diag$Index_3yr<-FALSE

	dat1diag$target_cancer_first<-1
	min_date_diagnosis<-unlist(lapply(1:nrow(dat1diag),FUN=function(x) 
			min(as.numeric(dat1diag[x,Names]),na.rm=TRUE)))
	dat1diag$first_diagnosis<-as.Date(min_date_diagnosis,origin="1970-01-01")
	dat<-rbind(dat1diag,dat1)
	return(dat)
}

	# diff_date_list

	# date_diag_list
	# class(ICD_list[751])
	# as.numeric(dat3$date_diagnosis2[100])
	# date_diag_list[100]
	

	# ICD_list[753]

	# 	if(all(is.na(ICD9)))
	# 	{
	# 		ICD[[i]]<-ICD10
	# 	}else{

	# 		stop("check icd9")
	# 	}

	# 	!is.na(ICD10)
	# 	dat3[43,c(Names[1:2],Names2[1:2],Names3[1:2])]
	# 	ICD9<-ICD9[!is.na(ICD9)]

	# 	Names_keep<-Names[1:length(Names) != Instance_min_date[[i]]] #exclude the instance of f.40005 with the earliest date of diagnosis for any cancer
	# 	min_date_diagnosis2_dat3[[i]]<-min(as.numeric(dat3[i,Names_keep]),na.rm=TRUE
	# 		)

	# 	# which instance contains the second earliest diagnosis date?
	# 	Instance_min_date2[[i]]<-which(as.numeric(dat3[i,Names]) == min_date_diagnosis2_dat3[[i]])
	# 	Names2_keep<-Names2[Instance_min_date2[[i]]]
	# 	if(sum(grep(icd10,dat3[i,Names2_keep]))

	# }
	

	# dat3$second_diagnosis<-as.Date(unlist(min_date_diagnosis2_dat3),origin="1970-01-01")

	# dat3[1:5,c("first_diagnosis","second_diagnosis","date_diagnosis2")]


	# Behaviour2<-unlist(lapply(1:nrow(dat),FUN=function(x)
	# 	unlist(Behaviour[x])[Instance_min_date[x]]))
	# dat$behaviour<-Behaviour2

	#create index to identify cases where lung cancer was first diagnosis. What is the minimum difference in time between the first diagnosis and any subsequent diagnosis (of non-target cancer). Our objective is e.g. to exclude individuals with multiple diagnoses and target cancer was the first diagnsosi but another cancer occured within 3/6 months of the first. The suscpicion here is that the target cancer is actually a metastasis from another site. e.g. individual presents with lung symptoms and diagnosis with lung cancer but this is secondary to another cancer that is only diagnosed later. 

	# 3.Exclude individuals with multiple diagnoses if LC is first diagnosis but the other cancer occurs within 3 months of the first (6 months?)
# 4.Exclude individuals with multiple diagnoses if LC is not first diagnosis and occurs <5 (3?) years after the first diagnosis



	# Test_list<-NULL
	# for(i in occurences){
	# 	Test<-dat[which(dat$N_diagnoses==i)[1],c("geneticID",Names,Names2,Names4,Names5,Names6)]
	# 	Pos<-unlist(lapply(1:ncol(Test),FUN=function(i) !is.na(Test[,i])))
	# 	Test[,Pos]
	# 	Test_list[[i]]<-Test[,Pos]
	# }
	# dat2<-do.call(plyr::rbind.fill,Test_list)
# 	return(dat2)
# }