
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UkbCancerMortality

<!-- badges: start -->

<!-- badges: end -->

The goal of UkbCancerMortality is to format UK Biobank cancer mortality
data so that it can be used in survival analyses. The package provides
functions for estimating the number of deaths due to the specific cancer
of interest and other causes, date of diagnosis, date of death,
difference between latest and earlies dates of diagnosis for individuals
with more than one diagnosis date, length of followup and survival time.

## Installation

You can install UkbCancerMortality from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("mightyphil2000/UkbCancerMortality")
```

## Example

In this example I load a dataset containing only lung cancer cases. I
then format the dataset to contain: to This is a basic example which
shows you how to solve a common problem:

``` r
library(UkbCancerMortality)
#load("../UKBB_cancer_mortality_incl_relateds_non_white_british.Rdata")
#lun<-bd[which(bd$overall_lung_cancer == 2),]
#lun2<-format_lung_cancer(dat=lun)
```
