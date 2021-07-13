# R function to extract the cancers from the cancer registry data
# Tom Dudding
# 25 Feb 2019

#* ICD10/9 code of cancer
#* Date of diagnosis
#* Histology of cancer
#* Behaviour of cancer (benign/malignant/metastatic)
#* Age at cancer diagnosis
#
#The function also removes duplicate entries (identified by date of diagnosis) **NB there may still be diagnoses that have very similar dates (these should be searched for manually if doing longitudinal analysis).**
#
#If a person has had more than one diagnosis of the cancer in question, this function creates a column for each separate diagnosis, this can be used to identify recurrence.
#
#By changing the option 'other' to FALSE the function will identify all cancers in the dataset **other than those of interest** which can be useful if wanting to select controls who are cancer free.

# Options
# dat - UK Biobank dataframe - the variable names should ideally match those provided by UK Biobank, NB the id column must be projectID
# cancerCode - vector of ICD codes of interest, for ICD10 this is a vector of class character, for ICD9 this is a vector of class integer.
# sitename - the label that will be given to the columns outputted. For a single diagnosis the column will be <sitename>1, subsequent diagnoses will be <sitename>2, <sitename>3 etc.
# other (default=F) - if this is set to true the function will identify all cancers that are not of interest (i.e. all except those in ICD10code) 
# cancer_col_start - the first column including the cancer data (the default corresponds to the default name of the ICD10 data)
# cancer_col_fin - the final colimn including the cancer data.
# date_col_start - first column of the cancer date variable (the same for both ICD10 and ICD9 so should not need to be changed from default)
# date_col_fin - last column of the cancer date variable
# hist_col_start - first column of the cancer histology variable (the same for both ICD10 and ICD9 so should not need to be changed from default)
# hist_col_fin - last column of the cancer histology variable
# behav_col_start - first column of the cancer behaviour variable (the same for both ICD10 and ICD9 so should not need to be changed from default)
# behav_col_fin - last column of the cancer behaviour variable
# age_col_start - first column of the cancer age variable (the same for both ICD10 and ICD9 so should not need to be changed from default)
# age_col_fin - last column of the cancer age variable

# To extract all cancer diagnoses (any ICD code) then select other=T and cancerCode=NA

# Output
# The function outputs the whole UKB dataset with the addition of a new column (labelled <sitename>#) for each incident cancer. In each column, for each person with a cancer there is a 
# character variable providing the information about the diagnosis in the following format:
# ICD_code/date_of_diagnosis/histology_code/behaviour_code/age_at_diagnosis
# There will be one column for each cancer diagnosis and these will be in order of earliest cancer first (e.g.sitename1 will be diagnosed before sitename2), duplicates have been removed 
# but if looking at multiple diagnoses it is still important to check the diagnoses and dates as some are very close together (likely referring to the same diagnosis) - this should be done manually.

UKBcancerFunc<-function(dat,cancerCode, sitename, other=F, 
                        cancer_col="f.40006.", 
                        date_col ="f.40005.", 
                        hist_col ="f.40011.",
                        behav_col = "f.40012.",
                        age_col = "f.40008.") {
 

  if (other ==T) {
    #generate list of ICDcodes not inlcuding those of interest
    allICD<-unlist(dplyr::select(dat,starts_with(cancer_col)), use.names = F)
    allICD<-allICD[!is.na(allICD)]
    allICD<-allICD[!duplicated(allICD)]
    cancerCode<-allICD[!allICD %in% cancerCode]
  }
  
  #idenitfy all participants with specific ICD10 codes
  cancers<-dplyr::select(dat,projectID,starts_with(cancer_col))
  cancerDate<-dplyr::select(dat,projectID,starts_with(date_col))
  cancerHist<-dplyr::select(dat,projectID,starts_with(hist_col))
  cancerBehav<-dplyr::select(dat,projectID,starts_with(behav_col))
  cancerAge<-dplyr::select(dat,projectID,starts_with(age_col))
  #rownames(cancers)<-rownames(dat)
  #rownames(cancerDate)<-rownames(dat)
  cancers2<-cancers[apply(cancers, 1, function(r) any(r %in% cancerCode)),]
  if (length(cancers2$projectID)<1) {
    message("No cancer of this type")
    return(dat)
    } else {
    
  cancerDate2<-cancerDate[apply(cancers, 1, function(r) any(r %in% cancerCode)),]
  cancerHist2<-cancerHist[apply(cancers, 1, function(r) any(r %in% cancerCode)),]
  cancerBehav2<-cancerBehav[apply(cancers, 1, function(r) any(r %in% cancerCode)),]
  cancerAge2<-cancerAge[apply(cancers, 1, function(r) any(r %in% cancerCode)),]
  
  #this loop takes the cancer codes and dates and for each participant generates a list of dates that they were diagnosed with the cancer of interest
  cancerList<-as.list(NULL)
  for (i in 1:length(cancers2$projectID)) {
    person1C<-unlist(cancers2[i,-1])
    person1D<-as.Date(unlist(cancerDate2[i,-1]),origin = "1970-01-01") #also converts to r internal dates (number of days since 1970-01-01)
    person1H<-unlist(cancerHist2[i,-1])
    person1B<-unlist(cancerBehav2[i,-1])
    person1A<-unlist(cancerAge2[i,-1])
    #set to missing non interest cancer codes
    person1C<-ifelse(person1C %in% cancerCode,person1C,NA) #keeps the NAs
    #person1C<-person1C[person1C %in% oral] #remves the NAs
    #remove dates/histology/behaviour/age for cancer codes not interested in
    person1D<-as.Date(ifelse(!is.na(person1C),person1D,NA),origin = "1970-01-01")
    person1D<-as.Date(person1D[!is.na(person1C)],origin = "1970-01-01")
    person1H<-ifelse(!is.na(person1C),person1H,NA)
    person1H<-person1H[!is.na(person1C)]
    person1B<-ifelse(!is.na(person1C),person1B,NA)
    person1B<-person1B[!is.na(person1C)]
    person1A<-ifelse(!is.na(person1C),person1A,NA)
    person1A<-person1A[!is.na(person1C)]
    person1C<-person1C[!is.na(person1C)]
    
    
    #remove duplicated cancers
    if (length(person1D)>1) {
      person1C<-person1C[!duplicated(person1D)]
      person1H<-person1H[!duplicated(person1D)]
      person1B<-person1B[!duplicated(person1D)]
      person1A<-person1A[!duplicated(person1D)]
      person1D<-as.Date(person1D[!duplicated(person1D)],origin = "1970-01-01")
      #sort by earliest diagnosis first
      iD <- order(person1D)
      person1D<-as.Date(person1D[iD],origin = "1970-01-01")
      person1C<-person1C[iD]
      person1H<-person1H[iD]
      person1B<-person1B[iD]
      person1A<-person1A[iD]
      
    } else {
      person1D<-as.Date(person1D,origin = "1970-01-01")
    }
    cancerList[[i]]<-paste(person1C,person1D,person1H,person1B,person1A, sep  = "/")
    #cancerList[[i]]<-c(person1C,as.Date(person1D,origin = "1970-01-01"),person1H,person1B,person1A)
  }
  
  #create a dataframe of cancer diagnosis dates
  
  cancerDF<-plyr::ldply(cancerList, rbind)
  maxC<-max(unlist(lapply(cancerList,length)))
  #cancerDF[,1:maxC]<-as.character(cancerDF[,1:maxC])
  colnames(cancerDF)<-paste0(sitename,1:maxC)

  
  cancerDF$projectID<-cancers2$projectID
  
  #merge cancer data back into fulldataset
  fulldat<-merge(dat,cancerDF, by="projectID", all.x = T)
  return(fulldat)
    }
}