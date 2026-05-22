#' Differential expression analysis using Schaum et al. (2020) Kidney data.
#'
#'  Expression data for mouse kidney samples across time
#'
#' @format A DESeq2 object containing results of differential expression analysis:
#' \describe{
#'   \item{assays}{List of matrices with raw and normalized count data}
#'   \item{rowRanges}{Genomic ranges for each gene}
#'   \item{colData}{Data frame with sample metadata}
#'   \item{design}{Formula representing the experimental design}
#'   \item{dispersions}{Estimated dispersion values}
#'   \item{results}{Data frame with differential expression results}
#' }
#' @source \url{doi:10.1038/s41586-020-2499-y}
"kidney"

#' Sample meta data for a longitudinal mouse kidney dataset
#'
#'  Sample meta data for the kidney data
#'
#' @format A data frame  with 30 rows and 3 variables:
#' \describe{
#'   \item{tissue}{name of organ the sample comes from}
#'   \item{age}{age of the mouse in months when sampled}
#'   \item{sex}{sex of mouse from which sample was taken}
#'   \item{age_sex}{combination of age and sex categorical variable}
#'   \item{tissue_age}{combination of tissue and age categorical variable}
#'   \item{condition}{combination of age, sex, and tissue categorical variable}
#'   \item{sample}{unique identifier for sample from original study}
#' }
#' @source \url{doi:10.1038/s41586-020-2499-y}
"kidney_metadata"

#' Differential expression analysis using Schaum et al. (2020) brain data.
#'
#'  Expression data for mouse brain samples across time
#'
#' @format A DESeq2 object containing results of differential expression analysis:
#' \describe{
#'   \item{assays}{List of matrices with raw and normalized count data}
#'   \item{rowRanges}{Genomic ranges for each gene}
#'   \item{colData}{Data frame with sample metadata}
#'   \item{design}{Formula representing the experimental design}
#'   \item{dispersions}{Estimated dispersion values}
#'   \item{results}{Data frame with differential expression results}
#' }
#' @source \url{doi:10.1038/s41586-020-2499-y}
"brain"

#' Sample metadata for a longitudinal Schaum et al. 2020 mouse brain dataset
#'
#'  Sample metadata for the Schaum et al. 2020 brain data
#'
#' @format A data frame  with 30 rows and 3 variables:
#' \describe{
#'   \item{tissue}{name of organ the sample comes from}
#'   \item{age}{age of the mouse in months when sampled}
#'   \item{sex}{sex of mouse from which sample was taken}
#'   \item{age_sex}{combination of age and sex categorical variable}
#'   \item{tissue_age}{combination of tissue and age categorical variable}
#'   \item{condition}{combination of age, sex, and tissue categorical variable}
#'   \item{sample}{unique identifier for sample from original study}
#' }
#' @source \url{doi:10.1038/s41586-020-2499-y}
"brain_metadata"



#' Differential expression analysis using Schaum et al. (2020) heart data.
#'
#'  Expression data for mouse heart samples across time
#'
#' @format A DESeq2 object containing results of differential expression analysis:
#' \describe{
#'   \item{assays}{List of matrices with raw and normalized count data}
#'   \item{rowRanges}{Genomic ranges for each gene}
#'   \item{colData}{Data frame with sample metadata}
#'   \item{design}{Formula representing the experimental design}
#'   \item{dispersions}{Estimated dispersion values}
#'   \item{results}{Data frame with differential expression results}
#' }
#' @source \url{doi:10.1038/s41586-020-2499-y}
"heart"

#' Sample metadata for a longitudinal Schaum et al. 2020 mouse heart dataset
#'
#'  Sample metadata for the Schaum et al. 2020 heart data
#'
#' @format A data frame  with 30 rows and 3 variables:
#' \describe{
#'   \item{tissue}{name of organ the sample comes from}
#'   \item{age}{age of the mouse in months when sampled}
#'   \item{sex}{sex of mouse from which sample was taken}
#'   \item{age_sex}{combination of age and sex categorical variable}
#'   \item{tissue_age}{combination of tissue and age categorical variable}
#'   \item{condition}{combination of age, sex, and tissue categorical variable}
#'   \item{sample}{unique identifier for sample from original study}
#' }
#' @source \url{doi:10.1038/s41586-020-2499-y}
"heart_metadata"

