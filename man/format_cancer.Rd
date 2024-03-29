% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_cancer.R
\name{format_cancer}
\alias{format_cancer}
\title{Format cancer mortality data}
\usage{
format_cancer(
  dat = NULL,
  censor_date = "2018-02-06",
  cancer_name = NULL,
  icd10 = NULL
)
}
\arguments{
\item{dat}{the cancer dataset}

\item{censor_date}{the date up until which the death registry data in UK Biobank is assumed to be complete. We assume that everyone not in the death registry at this date was still alive on this date. The default assumes a censoring date of "2018-02-06", based on the maximum value for the date of diagnosis f.40005.0.0 field in dataset 37205 on the rdsf /projects/MRC-IEU/research/data/ukbiobank/phenotypic/applications/15825/released/2019-05-02/data/derived/formats/r}

\item{cancer_name}{the name of the field/column in dat corresponding to case-control status for the cancer of interest. The field must take on values of 1 (control) or 2 (cases).}

\item{icd10}{the icd10 code for the cancer of interest. Icd codes are hierarchically structured. Therefore, it is not necessary to include all the lowest level icd codes for the cancer of interest. Only the highest level code should be specified. For example, breast cancer corresponds to icd10 codes C50.0-C50.9. Therefore, to create a breast cancer dataset, you would only need to set icd10 to "C50" (ie icd10="C50"). It is only necessary to include multiple icd codes when the cancer is defined by more than one high level code. E.g. colorectal cancer is defined as icd10 codes C18, C19 and C20, in which case you would specify icd10=c("C18","C19","C20").}
}
\value{
data frame
}
\description{
Format cancer mortality dataset so that it can be used in survival analyses. The
function estimates the number of deaths due to the target cancer of interest and other causes, date of diagnosis, date of death, difference between latest and earlies dates of diagnosis for individuals with more than one diagnosis date, length of followup and survival time. There are two instances of primary cause of death in UK Biobank (f.40001.0.0 and f.40001.1.0). If the second instance is different to the first instance (ie an individual has more than one primary cause of death), the function checks to see if either of these instances match the indicated cancer. If there is a match, then the primary cause of death is set to that cancer. E.g. a case has prostate cancer (C61) and Pulmonary edema (J81) listed as primary causes of death. In this case, we would set the primary cause of death to prostate cancer. This is a rare occurence. I've only observed this example once for prostate cancer and 0 times for lung cancer, breast cancer and overall cancer. If a case has two primary causes of death and neither one matches the indicated cancer, then the function stops and returns a warning message.
}
