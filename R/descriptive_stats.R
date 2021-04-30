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