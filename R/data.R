# ===================================================================
# DATA DOCUMENTATION SCRIPT
# This file contains the documentation for all datasets in the package.
# ===================================================================

#' Colon and Rectal Cancer Clinical and Virome Data
#'
#' A dataset containing cleaned key clinical, survival, and virome abundance data
#' for patients with colon and rectal adenocarcinoma from the TCGA PanCancer Atlas study.
#'
#' @format A data frame with patient data in rows and variables in columns:
#' \describe{
#'   \item{sample_id}{Unique TCGA patient identifier.}
#'   \item{time_months}{Overall survival time in months, censored at 60 months.}
#'   \item{status}{Overall survival status: 1 for deceased, 0 for alive.}
#'   \item{age}{Age of the patient at diagnosis.}
#'   \item{sex}{Sex of the patient.}
#'   \item{race}{Race of the patient.}
#'   \item{ethnicity}{Ethnicity of the patient.}
#'   \item{stage}{AJCC pathologic tumor stage.}
#'   \item{...}{The remaining 138 columns correspond to different viral taxa,
#'     with each column containing the normalised abundance for that virus.}
#' }
#' @source The data was sourced, cleaned, and merged from the supplementary
#'   materials of Smyth et al. (2024) \doi{10.1002/cam4.70434}, originating
#'   from the cBioPortal for Cancer Genomics ('TCGA, PanCancer Atlas' study).
"crc_virome"

#' Simulated Rapeseed Germination and Growth Data
#'
#' A dataset simulating the germination and early growth of rapeseed
#' (Brassica napus L.) across a range of temperatures. The simulation is based
#' on parameters published by Haj Sghaier et al. (2022).
#'
#' @format A data frame with 1040 rows and 5 variables:
#' \describe{
#'   \item{temperature}{The experimental temperature in degrees Celsius.}
#'   \item{replicate}{The experimental replicate number.}
#'   \item{time}{The day of measurement for each seed.}
#'   \item{growth}{A simulated measure of seedling growth.}
#'   \item{germinated}{The germination status: 1 for germinated, 0 for not germinated.}
#' }
#' @source This dataset was computationally simulated based on parameters derived from
#'   Haj Sghaier et al. (2022) \doi{10.3390/plants11212819}.
"germination"
