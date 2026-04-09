# ===================================================================
# DATA DOCUMENTATION SCRIPT
# Documentation for all package datasets.
# ===================================================================

#' Colon and Rectal Cancer Clinical and Virome Data
#'
#' Cleaned clinical, survival, and virome abundance data for patients
#' with colon and rectal adenocarcinoma (TCGA PanCancer Atlas).
#'
#' @srrstats {G5.1} Exported dataset provided for examples and testing.
#'
#' @format A data frame with patient data in rows and variables in columns:
#'   \describe{
#'     \item{sample_id}{Unique TCGA patient identifier (character).}
#'     \item{time_months}{Overall survival (mos), censored 60 mos. (numeric).}
#'     \item{status}{Overall survival: 1 = deceased, 0 = alive (integer).}
#'     \item{age}{Age at diagnosis (numeric).}
#'     \item{sex}{Sex of the patient (character).}
#'     \item{race}{Race of the patient (character).}
#'     \item{ethnicity}{Ethnicity of the patient (character).}
#'     \item{stage}{AJCC pathologic tumour stage (character).}
#'     \item{Alfamovirus}{Normalised abundance (numeric).}
#'     \item{Allexivirus}{Normalised abundance (numeric).}
#'     \item{Alphacoronavirus}{Normalised abundance (numeric).}
#'     \item{Alphaentomopoxvirus}{Normalised abundance (numeric).}
#'     \item{Alphafusellovirus}{Normalised abundance (numeric).}
#'     \item{Alphapapillomavirus}{Normalised abundance (numeric).}
#'     \item{Alphapartitivirus}{Normalised abundance (numeric).}
#'     \item{Alpharetrovirus}{Normalised abundance (numeric).}
#'     \item{Alphatorquevirus}{Normalised abundance (numeric).}
#'     \item{Alphavirus}{Normalised abundance (numeric).}
#'     \item{Ambidensovirus}{Normalised abundance (numeric).}
#'     \item{Ampelovirus}{Normalised abundance (numeric).}
#'     \item{Andromedalikevirus}{Normalised abundance (numeric).}
#'     \item{Aparavirus}{Normalised abundance (numeric).}
#'     \item{Aquamavirus}{Normalised abundance (numeric).}
#'     \item{Arenavirus}{Normalised abundance (numeric).}
#'     \item{Avastrovirus}{Normalised abundance (numeric).}
#'     \item{Bacilladnavirus}{Normalised abundance (numeric).}
#'     \item{Bafinivirus}{Normalised abundance (numeric).}
#'     \item{Barnyardlikevirus}{Normalised abundance (numeric).}
#'     \item{Batrachovirus}{Normalised abundance (numeric).}
#'     \item{Bcep22likevirus}{Normalised abundance (numeric).}
#'     \item{Betacoronavirus}{Normalised abundance (numeric).}
#'     \item{Betaentomopoxvirus}{Normalised abundance (numeric).}
#'     \item{Betapapillomavirus}{Normalised abundance (numeric).}
#'     \item{Betapartitivirus}{Normalised abundance (numeric).}
#'     \item{Betaretrovirus}{Normalised abundance (numeric).}
#'     \item{Betatorquevirus}{Normalised abundance (numeric).}
#'     \item{Bicaudavirus}{Normalised abundance (numeric).}
#'     \item{Bracovirus}{Normalised abundance (numeric).}
#'     \item{Bromovirus}{Normalised abundance (numeric).}
#'     \item{Cafeteriavirus}{Normalised abundance (numeric).}
#'     \item{Capripoxvirus}{Normalised abundance (numeric).}
#'     \item{Carlavirus}{Normalised abundance (numeric).}
#'     \item{Caulimovirus}{Normalised abundance (numeric).}
#'     \item{Cavemovirus}{Normalised abundance (numeric).}
#'     \item{Cervidpoxvirus}{Normalised abundance (numeric).}
#'     \item{Chlorovirus}{Normalised abundance (numeric).}
#'     \item{Closterovirus}{Normalised abundance (numeric).}
#'     \item{Comovirus}{Normalised abundance (numeric).}
#'     \item{Cp220likevirus}{Normalised abundance (numeric).}
#'     \item{Crinivirus}{Normalised abundance (numeric).}
#'     \item{Cripavirus}{Normalised abundance (numeric).}
#'     \item{Cucumovirus}{Normalised abundance (numeric).}
#'     \item{Cytomegalovirus}{Normalised abundance (numeric).}
#'     \item{Deltabaculovirus}{Normalised abundance (numeric).}
#'     \item{Dicipivirus}{Normalised abundance (numeric).}
#'     \item{Dyopipapillomavirus}{Normalised abundance (numeric).}
#'     \item{Emaravirus}{Normalised abundance (numeric).}
#'     \item{Endornavirus}{Normalised abundance (numeric).}
#'     \item{Enterovirus}{Normalised abundance (numeric).}
#'     \item{Fabavirus}{Normalised abundance (numeric).}
#'     \item{Flavivirus}{Normalised abundance (numeric).}
#'     \item{Furovirus}{Normalised abundance (numeric).}
#'     \item{Gammacoronavirus}{Normalised abundance (numeric).}
#'     \item{Gammapapillomavirus}{Normalised abundance (numeric).}
#'     \item{Gammaretrovirus}{Normalised abundance (numeric).}
#'     \item{Gammatorquevirus}{Normalised abundance (numeric).}
#'     \item{Hepacivirus}{Normalised abundance (numeric).}
#'     \item{Hepatovirus}{Normalised abundance (numeric).}
#'     \item{Higrevirus}{Normalised abundance (numeric).}
#'     \item{Hordeivirus}{Normalised abundance (numeric).}
#'     \item{Hypovirus}{Normalised abundance (numeric).}
#'     \item{I3likevirus}{Normalised abundance (numeric).}
#'     \item{Ichnovirus}{Normalised abundance (numeric).}
#'     \item{Iflavirus}{Normalised abundance (numeric).}
#'     \item{Ilarvirus}{Normalised abundance (numeric).}
#'     \item{Influenzavirus_C}{Normalised abundance (numeric).}
#'     \item{Isavirus}{Normalised abundance (numeric).}
#'     \item{Kobuvirus}{Normalised abundance (numeric).}
#'     \item{L5likevirus}{Normalised abundance (numeric).}
#'     \item{Lambdalikevirus}{Normalised abundance (numeric).}
#'     \item{Lambdapapillomavirus}{Normalised abundance (numeric).}
#'     \item{Lentivirus}{Normalised abundance (numeric).}
#'     \item{Leporipoxvirus}{Normalised abundance (numeric).}
#'     \item{Lymphocryptovirus}{Normalised abundance (numeric).}
#'     \item{Mamastrovirus}{Normalised abundance (numeric).}
#'     \item{Mardivirus}{Normalised abundance (numeric).}
#'     \item{Mastadenovirus}{Normalised abundance (numeric).}
#'     \item{Microvirus}{Normalised abundance (numeric).}
#'     \item{Mimivirus}{Normalised abundance (numeric).}
#'     \item{Molluscipoxvirus}{Normalised abundance (numeric).}
#'     \item{Muromegalovirus}{Normalised abundance (numeric).}
#'     \item{Muscavirus}{Normalised abundance (numeric).}
#'     \item{N4likevirus}{Normalised abundance (numeric).}
#'     \item{Narnavirus}{Normalised abundance (numeric).}
#'     \item{Nepovirus}{Normalised abundance (numeric).}
#'     \item{Nyavirus}{Normalised abundance (numeric).}
#'     \item{Omegapapillomavirus}{Normalised abundance (numeric).}
#'     \item{Orthobunyavirus}{Normalised abundance (numeric).}
#'     \item{Orthohepadnavirus}{Normalised abundance (numeric).}
#'     \item{Orthopoxvirus}{Normalised abundance (numeric).}
#'     \item{Ostreavirus}{Normalised abundance (numeric).}
#'     \item{Parapoxvirus}{Normalised abundance (numeric).}
#'     \item{Pecluvirus}{Normalised abundance (numeric).}
#'     \item{Pestivirus}{Normalised abundance (numeric).}
#'     \item{Phi29likevirus}{Normalised abundance (numeric).}
#'     \item{Phikmvlikevirus}{Normalised abundance (numeric).}
#'     \item{Phikzlikevirus}{Normalised abundance (numeric).}
#'     \item{Piscihepevirus}{Normalised abundance (numeric).}
#'     \item{Pithovirus}{Normalised abundance (numeric).}
#'     \item{Polyomavirus}{Normalised abundance (numeric).}
#'     \item{Pomovirus}{Normalised abundance (numeric).}
#'     \item{Potexvirus}{Normalised abundance (numeric).}
#'     \item{Potyvirus}{Normalised abundance (numeric).}
#'     \item{Prasinovirus}{Normalised abundance (numeric).}
#'     \item{Proboscivirus}{Normalised abundance (numeric).}
#'     \item{Protoparvovirus}{Normalised abundance (numeric).}
#'     \item{Prymnesiovirus}{Normalised abundance (numeric).}
#'     \item{Ranavirus}{Normalised abundance (numeric).}
#'     \item{Rhadinovirus}{Normalised abundance (numeric).}
#'     \item{Rubulavirus}{Normalised abundance (numeric).}
#'     \item{Salivirus}{Normalised abundance (numeric).}
#'     \item{Sapelovirus}{Normalised abundance (numeric).}
#'     \item{Sclerodarnavirus}{Normalised abundance (numeric).}
#'     \item{Senecavirus}{Normalised abundance (numeric).}
#'     \item{Sfi1unalikevirus}{Normalised abundance (numeric).}
#'     \item{Sfi21dtunalikevirus}{Normalised abundance (numeric).}
#'     \item{Sicinivirus}{Normalised abundance (numeric).}
#'     \item{Simplexvirus}{Normalised abundance (numeric).}
#'     \item{Skunalikevirus}{Normalised abundance (numeric).}
#'     \item{Sp6likevirus}{Normalised abundance (numeric).}
#'     \item{Spo1virus}{Normalised abundance (numeric).}
#'     \item{Spounalikevirus}{Normalised abundance (numeric).}
#'     \item{T4likevirus}{Normalised abundance (numeric).}
#'     \item{T5likevirus}{Normalised abundance (numeric).}
#'     \item{Taupapillomavirus}{Normalised abundance (numeric).}
#'     \item{Tenuivirus}{Normalised abundance (numeric).}
#'     \item{Tobamovirus}{Normalised abundance (numeric).}
#'     \item{Totivirus}{Normalised abundance (numeric).}
#'     \item{Trichovirus}{Normalised abundance (numeric).}
#'     \item{Tritimovirus}{Normalised abundance (numeric).}
#'     \item{Tunalikevirus}{Normalised abundance (numeric).}
#'     \item{Tymovirus}{Normalised abundance (numeric).}
#'     \item{Vesivirus}{Normalised abundance (numeric).}
#'     \item{Waikavirus}{Normalised abundance (numeric).}
#'     \item{Whispovirus}{Normalised abundance (numeric).}
#'     \item{Yatapoxvirus}{Normalised abundance (numeric).}
#'   }
#'
#' @source Sourced, cleaned, and merged from Smyth et al. (2024)
#'   [\doi{10.1002/cam4.70434}], originating from cBioPortal
#'   (TCGA, PanCancer Atlas).
#'   \url{https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018}
#'
#' @references
#' Smyth, J., Godet, J., Choudhary, A., Das, A., Gkoutos, G.
#' V., & Acharjee, A. (2024). Microbiome-Based
#' Colon Cancer Patient Stratification and Survival Analysis. *Cancer
#' Medicine*,13(22), e70434. https://doi.org/10.1002/cam4.70434
#'
#' @keywords datasets
#' @name crc_virome
NULL

#' Simulated Rapeseed Germination and Growth Data
#'
#' Simulated germination and early growth of rapeseed (Brassica napus L.)
#' across a range of temperatures. Based on parameters from
#' Haj Sghaier et al. (2022).
#'
#' @srrstats {G5.1} Exported dataset provided for examples and testing.
#'
#' @format A data frame with 1040 rows and 5 variables:
#' \describe{
#'   \item{temperature}{Experimental temperature (Celsius) (numeric).}
#'   \item{replicate}{Experimental replicate number (integer).}
#'   \item{time}{Day of measurement (integer).}
#'   \item{growth}{Simulated seedling growth (numeric).}
#'   \item{germinated}{Germination status: 1 = germinated,
#'      0 = not germinated (integer).}
#' }
#'
#' @source Computationally simulated based on parameters derived from
#'   Haj Sghaier et al. (2022) \doi{10.3390/plants11212819}.
#'
#' @references
#' Haj Sghaier, A., Tarnawa, Á., Khaeim, H., Kovács, G. P.,
#' Gyuricza, C., & Kende, Z. (2022). The Effects of Temperature and
#' Water on the Seed Germination and Seedling Development of Rapeseed
#' (Brassica napus L.). *Plants*, 11(21), 2819.
#' \doi{10.3390/plants11212819}
#'
#' @keywords datasets
#' @name germination
NULL
