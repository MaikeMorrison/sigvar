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
#'   \item{mSBS1, mSBS5, ..., mSBS_N3}{Number of mutations for each signature}
#'   \item{Tissue}{Tissue of origin}
#'   \item{DIAGNOSIS}{Type of cancer}
#'   \item{dose}{Treatment dose}
#'   \item{chemical}{Known or suspected carcinogen chemical administered}
#' }
#' @source <https://github.com/team113sanger/mouse-mutatation-signatures/blob/6a00d910df40d178c373ac4d57849918e9814951/figure1/mexposure.rds>
"mutsig_carcinogens_mice_SBS"




#' Bootstrap p-values comparing SVA results for carcinogen-exposed murine tumors to spontaneous tumors
#'
#' Data on single Base Substitution mutational signature attributions
#' from Riva et al. Nat Genet 52, 1189–1197 (2020). https://doi.org/10.1038/s41588-020-0692-4
#'
#' @format ## `mutsig_carcinogens_mice_bootstrap_p_vals`
#' A tibble with 29 rows and 7 columns:
#' \describe{
#'   \item{group}{Organ of tumor and exposure status}
#'   \item{chemical}{Exposure substance}
#'   \item{Tissue}{Tissue of origin}
#'   \item{group_2}{Tumor/exposure reference each tumor/exposure group was compared to (either "Spontaneous_liver"  or"Spontaneous_lung" )}
#'   \item{across_sample_heterogeneity}{P-values comparing values of across-sample heterogeneity between group and group_2}
#'   \item{mean_within_sample_diversity}{P-values comparing values of within-sample diversity between group and group_2}
#'   \item{pooled_diversity}{P-values comparing the pooled diversity between group and group_2}
#' }
"mutsig_carcinogens_mice_bootstrap_p_vals"




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

#' Color palette for mutational signatures analyzed in the study by Morrison et al.
#'
#' Vector of named colors for use in plotting the relative activities of mutational
#' signatures.
#'
#' @format ## `all_sig_pal`
#' A vector of colors for mutational signature activities
#' \describe{
#' Hex color codes for mutational signatures
#' }
"all_sig_pal"


#' Activity of SBS, DBS, and ID mutational signatures in a global sample of esophageal squamous cell carcinoma cases.
#'
#' Data from: Moody, S., Senkin, S., Islam, S.M.A. et al. Mutational signatures in esophageal squamous cell carcinoma from eight countries with varying incidence. Nat Genet 53, 1553–1563 (2021).
#'
#' @format ## `ESCC_sig_activity`
#' A dataframe with 552 rows and 46 columns:
#' \describe{
#' \item{Country}{Country where each sample was collected}
#' \item{Incidence_Level}{Incidence level (High or Low) of ESCC in the region where the sample was collected}
#' \item{Sample}{Unique sample ID}
#' \item{SBS1, ..., ID17}{Relative abundance of mutational signatures in each sample. These columns sum to 1 for each row.}
#' }
#' @source <https://doi.org/10.1038/s41588-021-00928-6>
"ESCC_sig_activity"


#' Cosine similarity of SBS, DBS, and ID mutational signatures found in a global sample of esophageal squamous cell carcinoma cases.
#'
#' Data from: Moody, S., Senkin, S., Islam, S.M.A. et al. Mutational signatures in esophageal squamous cell carcinoma from eight countries with varying incidence. Nat Genet 53, 1553–1563 (2021).
#'
#' @format ## `ESCC_sig_similarity`
#' A dataframe with 43 rows and 43 columns:
#' \describe{
#' \item{SBS1, ..., ID17}{Pairwise cosine similarities between mutational signatures.}
#' }
#' @source <https://doi.org/10.1038/s41588-021-00928-6>
"ESCC_sig_similarity"





#' Activity of SBS mutational signatures in samples of lung cancer in never smokers
#'
#' Data from: Zhang et al. (2021). Genomic and evolutionary classification of lung cancer in never smokers. Nature Genetics, 53(9), 1348–1359. https://doi.org/10.1038/s41588-021-00920-0
#'
#' @format ## `zhang_sig_activity `
#' A dataframe with 39 rows and 213 columns:
#' \describe{
#' \item{Subject}{Unique ID for each subject}
#' ...
#' \item{passive_smoking}{Subject smoking status ("Non-smoker" or "Passive smoker")}
#' ...
#' \item{SBS1, ..., SBS40}{Relative abundance of mutational signatures in each sample. These columns sum to 1 for each row.}
#' }
#' @source <https://doi.org/10.1038/s41588-021-00920-0>
"zhang_sig_activity "


#' Cosine similarity of SBS mutational signatures found in samples of lung cancer in never smokers
#'
#' Data from: Zhang et al. (2021). Genomic and evolutionary classification of lung cancer in never smokers. Nature Genetics, 53(9), 1348–1359. https://doi.org/10.1038/s41588-021-00920-0

#'
#' @format ## `zhang_sig_similarity`
#' A dataframe with 14 rows and 14 columns:
#' \describe{
#' \item{SBS1, ..., SBS40}{Pairwise cosine similarities between mutational signatures.}
#' }
#' @source <https://doi.org/10.1038/s41588-021-00920-0>
"zhang_sig_similarity"



#' Activity of Smoking, APOBEC, and Ageing mutational signatures in a sample of East Asian lung adenocarcinoma patients.
#'
#' Data from: Chen, J., Yang, H., Teo, A.S.M. et al. Genomic landscape of lung adenocarcinoma in East Asians. Nat Genet 52, 177–186 (2020).
#'
#' @format ## `smoker_sigs_chen`
#' A dataframe with 88 rows and 20 columns:
#' \describe{
#' \item{Patient.ID}{Unique identifier for each patient}
#' \item{Smoker}{Smoking status for each patient: "Smoker" or "Non-smoker"}
#' \item{Cohort, Stage, Age, Gender}{Patient metadata}
#' \item{Smoking}{Activity of smoking mutational signature}
#' \item{APOBEC}{Activity of APOBEC mutational signature}
#' \item{Ageing}{Activity of ageing mutational signature}
#' }
#' @source <https://doi.org/10.1038/s41588-019-0569-6>
"smoker_sigs_chen"


#' Cosine similarity of smoking, APOBEC, and ageing mutational signatures found in a sample of East Asian lung adenocarcinoma patients.
#'
#' Data from: Chen, J., Yang, H., Teo, A.S.M. et al. Genomic landscape of lung adenocarcinoma in East Asians. Nat Genet 52, 177–186 (2020).
#'
#' @format ## `ESCC_sig_similarity`
#' A dataframe with 3 rows and 3 columns:
#' \describe{
#' \item{Smoking, ..., Ageing}{Pairwise cosine similarities between mutational signatures.}
#' }
#' @source <https://doi.org/10.1038/s41588-019-0569-6>
"ESCC_sig_similarity"
