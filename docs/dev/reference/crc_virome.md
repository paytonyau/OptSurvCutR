# Colon and Rectal Cancer Clinical and Virome Data

Cleaned clinical, survival, and virome abundance data for patients with
colon and rectal adenocarcinoma (TCGA PanCancer Atlas).

## Format

A data frame with patient data in rows and variables in columns:

- sample_id:

  Unique TCGA patient identifier (character).

- time_months:

  Overall survival (months), censored 60 mos. (numeric).

- status:

  Overall survival status: 1 = deceased, 0 = alive (integer).

- age:

  Age at diagnosis (numeric).

- sex:

  Sex of the patient (character).

- race:

  Race of the patient (character).

- ethnicity:

  Ethnicity of the patient (character).

- stage:

  AJCC pathologic tumour stage (character).

- Alfamovirus:

  Normalised abundance (numeric).

- Allexivirus:

  Normalised abundance (numeric).

- Alphacoronavirus:

  Normalised abundance (numeric).

- Alphaentomopoxvirus:

  Normalised abundance (numeric).

- Alphafusellovirus:

  Normalised abundance (numeric).

- Alphapapillomavirus:

  Normalised abundance (numeric).

- Alphapartitivirus:

  Normalised abundance (numeric).

- Alpharetrovirus:

  Normalised abundance (numeric).

- Alphatorquevirus:

  Normalised abundance (numeric).

- Alphavirus:

  Normalised abundance (numeric).

- Ambidensovirus:

  Normalised abundance (numeric).

- Ampelovirus:

  Normalised abundance (numeric).

- Andromedalikevirus:

  Normalised abundance (numeric).

- Aparavirus:

  Normalised abundance (numeric).

- Aquamavirus:

  Normalised abundance (numeric).

- Arenavirus:

  Normalised abundance (numeric).

- Avastrovirus:

  Normalised abundance (numeric).

- Bacilladnavirus:

  Normalised abundance (numeric).

- Bafinivirus:

  Normalised abundance (numeric).

- Barnyardlikevirus:

  Normalised abundance (numeric).

- Batrachovirus:

  Normalised abundance (numeric).

- Bcep22likevirus:

  Normalised abundance (numeric).

- Betacoronavirus:

  Normalised abundance (numeric).

- Betaentomopoxvirus:

  Normalised abundance (numeric).

- Betapapillomavirus:

  Normalised abundance (numeric).

- Betapartitivirus:

  Normalised abundance (numeric).

- Betaretrovirus:

  Normalised abundance (numeric).

- Betatorquevirus:

  Normalised abundance (numeric).

- Bicaudavirus:

  Normalised abundance (numeric).

- Bracovirus:

  Normalised abundance (numeric).

- Bromovirus:

  Normalised abundance (numeric).

- Cafeteriavirus:

  Normalised abundance (numeric).

- Capripoxvirus:

  Normalised abundance (numeric).

- Carlavirus:

  Normalised abundance (numeric).

- Caulimovirus:

  Normalised abundance (numeric).

- Cavemovirus:

  Normalised abundance (numeric).

- Cervidpoxvirus:

  Normalised abundance (numeric).

- Chlorovirus:

  Normalised abundance (numeric).

- Closterovirus:

  Normalised abundance (numeric).

- Comovirus:

  Normalised abundance (numeric).

- Cp220likevirus:

  Normalised abundance (numeric).

- Crinivirus:

  Normalised abundance (numeric).

- Cripavirus:

  Normalised abundance (numeric).

- Cucumovirus:

  Normalised abundance (numeric).

- Cytomegalovirus:

  Normalised abundance (numeric).

- Deltabaculovirus:

  Normalised abundance (numeric).

- Dicipivirus:

  Normalised abundance (numeric).

- Dyopipapillomavirus:

  Normalised abundance (numeric).

- Emaravirus:

  Normalised abundance (numeric).

- Endornavirus:

  Normalised abundance (numeric).

- Enterovirus:

  Normalised abundance (numeric).

- Fabavirus:

  Normalised abundance (numeric).

- Flavivirus:

  Normalised abundance (numeric).

- Furovirus:

  Normalised abundance (numeric).

- Gammacoronavirus:

  Normalised abundance (numeric).

- Gammapapillomavirus:

  Normalised abundance (numeric).

- Gammaretrovirus:

  Normalised abundance (numeric).

- Gammatorquevirus:

  Normalised abundance (numeric).

- Hepacivirus:

  Normalised abundance (numeric).

- Hepatovirus:

  Normalised abundance (numeric).

- Higrevirus:

  Normalised abundance (numeric).

- Hordeivirus:

  Normalised abundance (numeric).

- Hypovirus:

  Normalised abundance (numeric).

- I3likevirus:

  Normalised abundance (numeric).

- Ichnovirus:

  Normalised abundance (numeric).

- Iflavirus:

  Normalised abundance (numeric).

- Ilarvirus:

  Normalised abundance (numeric).

- Influenzavirus_C:

  Normalised abundance (numeric).

- Isavirus:

  Normalised abundance (numeric).

- Kobuvirus:

  Normalised abundance (numeric).

- L5likevirus:

  Normalised abundance (numeric).

- Lambdalikevirus:

  Normalised abundance (numeric).

- Lambdapapillomavirus:

  Normalised abundance (numeric).

- Lentivirus:

  Normalised abundance (numeric).

- Leporipoxvirus:

  Normalised abundance (numeric).

- Lymphocryptovirus:

  Normalised abundance (numeric).

- Mamastrovirus:

  Normalised abundance (numeric).

- Mardivirus:

  Normalised abundance (numeric).

- Mastadenovirus:

  Normalised abundance (numeric).

- Microvirus:

  Normalised abundance (numeric).

- Mimivirus:

  Normalised abundance (numeric).

- Molluscipoxvirus:

  Normalised abundance (numeric).

- Muromegalovirus:

  Normalised abundance (numeric).

- Muscavirus:

  Normalised abundance (numeric).

- N4likevirus:

  Normalised abundance (numeric).

- Narnavirus:

  Normalised abundance (numeric).

- Nepovirus:

  Normalised abundance (numeric).

- Nyavirus:

  Normalised abundance (numeric).

- Omegapapillomavirus:

  Normalised abundance (numeric).

- Orthobunyavirus:

  Normalised abundance (numeric).

- Orthohepadnavirus:

  Normalised abundance (numeric).

- Orthopoxvirus:

  Normalised abundance (numeric).

- Ostreavirus:

  Normalised abundance (numeric).

- Parapoxvirus:

  Normalised abundance (numeric).

- Pecluvirus:

  Normalised abundance (numeric).

- Pestivirus:

  Normalised abundance (numeric).

- Phi29likevirus:

  Normalised abundance (numeric).

- Phikmvlikevirus:

  Normalised abundance (numeric).

- Phikzlikevirus:

  Normalised abundance (numeric).

- Piscihepevirus:

  Normalised abundance (numeric).

- Pithovirus:

  Normalised abundance (numeric).

- Polyomavirus:

  Normalised abundance (numeric).

- Pomovirus:

  Normalised abundance (numeric).

- Potexvirus:

  Normalised abundance (numeric).

- Potyvirus:

  Normalised abundance (numeric).

- Prasinovirus:

  Normalised abundance (numeric).

- Proboscivirus:

  Normalised abundance (numeric).

- Protoparvovirus:

  Normalised abundance (numeric).

- Prymnesiovirus:

  Normalised abundance (numeric).

- Ranavirus:

  Normalised abundance (numeric).

- Rhadinovirus:

  Normalised abundance (numeric).

- Rubulavirus:

  Normalised abundance (numeric).

- Salivirus:

  Normalised abundance (numeric).

- Sapelovirus:

  Normalised abundance (numeric).

- Sclerodarnavirus:

  Normalised abundance (numeric).

- Senecavirus:

  Normalised abundance (numeric).

- Sfi1unalikevirus:

  Normalised abundance (numeric).

- Sfi21dtunalikevirus:

  Normalised abundance (numeric).

- Sicinivirus:

  Normalised abundance (numeric).

- Simplexvirus:

  Normalised abundance (numeric).

- Skunalikevirus:

  Normalised abundance (numeric).

- Sp6likevirus:

  Normalised abundance (numeric).

- Spo1virus:

  Normalised abundance (numeric).

- Spounalikevirus:

  Normalised abundance (numeric).

- T4likevirus:

  Normalised abundance (numeric).

- T5likevirus:

  Normalised abundance (numeric).

- Taupapillomavirus:

  Normalised abundance (numeric).

- Tenuivirus:

  Normalised abundance (numeric).

- Tobamovirus:

  Normalised abundance (numeric).

- Totivirus:

  Normalised abundance (numeric).

- Trichovirus:

  Normalised abundance (numeric).

- Tritimovirus:

  Normalised abundance (numeric).

- Tunalikevirus:

  Normalised abundance (numeric).

- Tymovirus:

  Normalised abundance (numeric).

- Vesivirus:

  Normalised abundance (numeric).

- Waikavirus:

  Normalised abundance (numeric).

- Whispovirus:

  Normalised abundance (numeric).

- Yatapoxvirus:

  Normalised abundance (numeric).

## Source

Sourced, cleaned, and merged from Smyth et al. (2024)
\[[doi:10.1002/cam4.70434](https://doi.org/10.1002/cam4.70434) \],
originating from cBioPortal (TCGA, PanCancer Atlas).
<https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018>

## srrstats compliance

.

## References

Smyth, J., Godet, J., Choudhary, A., Das, A., Gkoutos, G. V., &
Acharjee, A. (2024). Microbiome-Based Colon Cancer Patient
Stratification and Survival Analysis. \*Cancer Medicine\*,13(22),
e70434. https://doi.org/10.1002/cam4.70434
