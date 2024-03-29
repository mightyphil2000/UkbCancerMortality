% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_lung_cancer.R
\name{format_lung_cancer}
\alias{format_lung_cancer}
\title{Format lung cancer mortality data}
\usage{
format_lung_cancer(dat = NULL, censor_date = "2018-02-06")
}
\arguments{
\item{dat}{the lung cancer dataset containing only individuals with a diagnosis of lung cancer}

\item{censor_date}{the date up until which the death registry data in UK Biobank is assumed to be complete. We assume that everyone not in the death registry at this date was still alive on this date. The default assumes a censoring date of "2018-02-06", based on the maximum value for the date of diagnosis f.40005.0.0 field in dataset 37205 on the rdsf /projects/MRC-IEU/research/data/ukbiobank/phenotypic/applications/15825/released/2019-05-02/data/derived/formats/r}
}
\value{
data frame
}
\description{
Format the lung cancer mortality dataset so that it can be used in survival analyses. The
function estimates the number of deaths due to lung cancer and other causes, date of diagnosis, date of death, difference between latest and earlies dates of diagnosis for individuals with more than one diagnosis date, length of followup and survival time.
}
