#' Single Base Substitution COSMIC mutational signatures for 2780 tumors
#'
#' Single Base Substitution mutational signature attributions to COSMIC signatures
#' from the Pan Cancer Whole-Genome consortium (PCAWG). Oct 19 2019 version
#'
#' @format ## `PCAWG_SigProfiler_COSMIC_SBS`
#' A tibble with 2,780 rows and 68 columns:
#' \describe{
#'   \item{Cancer Types}{Type of cancer}
#'   \item{Sample Names}{Identifier of the tumor}
#'   \item{SBS1, SBS2}{Number of mutations for each signature}
#'   ...
#' }
#' @source <https://dcc.icgc.org/releases/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples>
"PCAWG_SigProfiler_COSMIC_SBS"

#' Indel COSMIC mutational signatures for 2780 tumors
#'
#' Indel mutational signature attributions to COSMIC signatures
#' from the Pan Cancer Whole-Genome consortium (PCAWG). Oct 19 2019 version
#'
#' @format ## `PCAWG_SigProfiler_COSMIC_ID`
#' A tibble with 2,780 rows and 20 columns:
#' \describe{
#'   \item{Cancer Types}{Type of cancer}
#'   \item{Sample Names}{Identifier of the tumor}
#'   \item{ID1, ID2}{Number of mutations for each signature}
#'   ...
#' }
#' @source <https://dcc.icgc.org/releases/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples>
"PCAWG_SigProfiler_COSMIC_ID"

#' Double Base Substitution COSMIC mutational signatures for 2780 tumors
#'
#' Double Base Substitution mutational signature attributions to COSMIC signatures
#' from the Pan Cancer Whole-Genome consortium (PCAWG). Oct 19 2019 version
#'
#' @format ## `PCAWG_SigProfiler_COSMIC_DBS`
#' A tibble with 2,670 rows and 14 columns:
#' \describe{
#'   \item{Cancer Types}{Type of cancer}
#'   \item{Sample Names}{Identifier of the tumor}
#'   \item{DBS1, DBS2}{Number of mutations for each signature}
#'   ...
#' }
#' @source <https://dcc.icgc.org/releases/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples>
"PCAWG_SigProfiler_COSMIC_DBS"


#' Single Base Substitution mutational signatures for 181 mice exposed to known or suspected carcinogens
#'
#' Single Base Substitution mutational signature attributions
#' from Riva et al. Nat Genet 52, 1189–1197 (2020). https://doi.org/10.1038/s41588-020-0692-4
#'
#' @format ## `mutsig_carcinogens_mice_SBS`
#' A tibble with 181 rows and 46 columns:
#' \describe{
#'   \item{Sample}{Identifier of the tumor}
#'   \item{mSBS1, mSBS5}{Number of mutations for each signature}
#'   ...
#'   \item{Tissue}{Tissue of origin}
#'   \item{DIAGNOSIS}{Type of cancer}
#'   \item{dose}{Treatment dose}
#'   \item{chemical}{Known or suspected carcinogen chemical administered}
#' }
#' @source <https://github.com/team113sanger/mouse-mutatation-signatures/blob/6a00d910df40d178c373ac4d57849918e9814951/figure1/mexposure.rds>
"mutsig_carcinogens_mice_SBS"


#' Copy Number and Single Base Substitution mutational signatures for 120 malignant pleural mesothelioma from the MESOMICS cohort
#'
#' Copy Number and Single Base Substitution mutational signature attributions
#' from Mangiante et al. Nat Genet 55, 607–618 (2023). https://doi.org/10.1038/s41588-023-01321-1
#'
#' @format ## `MESOMICS_CN_SBS`
#' A tibble with 120 rows and 41 columns:
#' \describe{
#'   \item{ID_MESOMICS}{Identifier of the patient}
#'   \item{Sample}{Identifier of the tumor sample}
#'   \item{Sex}{Patient sex}
#'   \item{Age}{Patient age at diagnosis}
#'   \item{Smoking, Professional.Asbestos, ...}{Patient exposure status}
#'   \item{Survival.Censor, Survival.Time}{Survival status (dead or alive) and time in months}
#'   \item{Asbestos.Exposure.probability, Asbestos.Exposure.frequency, ...}{Quantitative asbestos exposure variables}
#'   \item{Type}{Tumor histological type (WHO)}
#'   \item{CNB}{Copy number burden}
#'   \item{CN1,...,CN19}{Copy number signature attributions}
#'   \item{SBS1,...,SBS40}{Single Base Substitution signature attributions}
#' }
#' @source <https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01321-1/MediaObjects/41588_2023_1321_MOESM4_ESM.xlsx>
"MESOMICS_CN_SBS"

#' Copy Number reference signatures for COSMIC3.1
#'
#' COSMIC 3.1 Copy Number mutational signatures
#' from Mangiante et al. Nat Genet 55, 607–618 (2023). https://doi.org/10.1038/s41588-023-01321-1
#'
#' @format ## `COSMIC3.1_CN`
#' A tibble with 48 rows and 8 columns:
#' \describe{
#'   \item{MutationsType}{Type of segment}
#'   \item{CN1,...,CN19}{Proportion of Copy number segments in each category}
#' }
#' @source <https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01321-1/MediaObjects/41588_2023_1321_MOESM4_ESM.xlsx>
"COSMIC3.1_CN"

#' Single Base Substitution COSMIC mutational signatures for 2780 tumors
#'
#' Single Base Substitution mutational signature attributions to COSMIC signatures
#' from the Pan Cancer Whole-Genome consortium (PCAWG). Oct 19 2019 version
#'
#' @format ## `sbs_palette`
#' A vector of colors for COSMIC SBS signatures
#' \describe{
#' RGB codes for COSMIC signatures
#' }
#' @source <https://dcc.icgc.org/releases/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples>
"sbs_palette"
