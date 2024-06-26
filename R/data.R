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

#' Single base substitution reference signatures for COSMIC3.0
#'
#' COSMIC 3.0 Single base subtitutions mutational signatures
#' from the COSMIC website
#'
#' @format ## `COSMIC3.0_SBS`
#' A tibble with 96 rows and 68 columns:
#' \describe{
#'   \item{Type}{Type of SBS}
#'   \item{SBS1,...,SBS85}{Proportion of SBS in each class for each signature}
#' }
#' @source <https://cancer.sanger.ac.uk/signatures/downloads/>
"COSMIC3.0_SBS"

#' Single base substitution reference signatures for COSMIC3.3.1
#'
#' COSMIC 3.3.1 Single base subtitutions mutational signatures
#' from the COSMIC website
#'
#' @format ## `COSMIC3.3.1_SBS`
#' A tibble with 96 rows and 80 columns:
#' \describe{
#'   \item{Type}{Type of SBS}
#'   \item{SBS1,...,SBS95}{Proportion of SBS in each class for each signature}
#' }
#' @source <https://cancer.sanger.ac.uk/signatures/downloads/>
"COSMIC3.3.1_SBS"

#' Indel reference signatures for COSMIC3.0
#'
#' COSMIC 3.0 Indel mutational signatures
#' from the COSMIC website
#'
#' @format ## `COSMIC3.0_ID`
#' A tibble with 83 rows and 18 columns:
#' \describe{
#'   \item{Type}{Type of ID}
#'   \item{ID1,...,ID17}{Proportion of ID in each class for each signature}
#' }
#' @source <https://cancer.sanger.ac.uk/signatures/downloads/>
"COSMIC3.0_ID"

#' Indel reference signatures for COSMIC3.3.1
#'
#' COSMIC 3.3 Indel mutational signatures
#' from the COSMIC website
#'
#' @format ## `COSMIC3.3_ID`
#' A tibble with 83 rows and 19 columns:
#' \describe{
#'   \item{Type}{Type of ID}
#'   \item{ID1,...,ID18}{Proportion of ID in each class for each signature}
#' }
#' @source <https://cancer.sanger.ac.uk/signatures/downloads/>
"COSMIC3.3_ID"

#' Double base reference signatures for COSMIC3.0
#'
#' COSMIC 3.3 Double base mutational signatures
#' from the COSMIC website
#'
#' @format ## `COSMIC3.3_DBS`
#' A tibble with 78 rows and 12 columns:
#' \describe{
#'   \item{Type}{Type of ID}
#'   \item{DBS1,...,DBS11}{Proportion of DBS in each class for each signature}
#' }
#' @source <https://cancer.sanger.ac.uk/signatures/downloads/>
"COSMIC3.3_DBS"

#' Mice Double base reference signatures for COSMIC3.0
#'
#' Mice COSMIC 3.3 Double base mutational signatures
#' from the COSMIC website
#'
#' @format ## `COSMIC3.3_DBS_mm10`
#' A tibble with 78 rows and 12 columns:
#' \describe{
#'   \item{Type}{Type of ID}
#'   \item{DBS1,...,DBS11}{Proportion of DBS in each class for each signature}
#' }
#' @source <https://cancer.sanger.ac.uk/signatures/downloads/>
"COSMIC3.3_DBS_mm10"

#' Mice Single base substitution reference signatures for COSMIC3.3.1
#'
#' Mice COSMIC 3.3.1 Single base subtitutions mutational signatures
#' from the COSMIC website
#'
#' @format ## `COSMIC3.3.1_SBS_mm10`
#' A tibble with 96 rows and 80 columns:
#' \describe{
#'   \item{Type}{Type of SBS}
#'   \item{SBS1,...,SBS95}{Proportion of SBS in each class for each signature}
#' }
#' @source <https://cancer.sanger.ac.uk/signatures/downloads/>
"COSMIC3.3.1_SBS_mm10"

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

#' Copy Number reference signatures for COSMIC3.3
#'
#' COSMIC 3.3 Copy Number mutational signatures
#' from the COSMIC website
#'
#' @format ## `COSMIC3.3_CN`
#' A tibble with 48 rows and 8 columns:
#' \describe{
#'   \item{MutationsType}{Type of segment}
#'   \item{CN1,...,CN19}{Proportion of Copy number segments in each category}
#' }
#' @source <https://cancer.sanger.ac.uk/signatures/downloads/>
"COSMIC3.3_CN"

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



#' Activity of SBS, DBS, and ID mutational signatures in a global sample of renal cell carcinoma cases.
#'
#' Data from: Senkin et al. (2023). Geographic variation of mutagenic exposures in kidney cancer genomes. MedRxiv, 14, 2023.06.20.23291538. https://doi.org/10.1101/2023.06.20.23291538
#'
#' @format ## `rcc_sig_activity`
#' A dataframe with 962 rows and 31 columns:
#' \describe{
#' \item{Sample}{Unique sample ID}
#' \item{Country}{Country where each sample was collected}
#' \item{ASR}{Age-standardized rate (ASR) of RCC incidence in each country}
#' \item{SBS1, ..., ID_C}{Relative abundance of mutational signatures in each sample. These columns sum to 1 for each row.}
#' }
#' @source <https://doi.org/10.1101/2023.06.20.23291538>
"rcc_sig_activity"


#' Activity of SBS, DBS, and ID mutational signatures in a global sample of renal cell carcinoma cases.
#'
#' Data from: Senkin et al. (2023). Geographic variation of mutagenic exposures in kidney cancer genomes. MedRxiv, 14, 2023.06.20.23291538. https://doi.org/10.1101/2023.06.20.23291538
#'
#' @format ## `rcc_sim`
#' A dataframe with 28 rows and 28 columns:
#' \describe{
#' \item{SBS1, ..., ID_C}{Cosine similarity between mutational signatures.}
#' }
#' @source <https://doi.org/10.1101/2023.06.20.23291538>
"rcc_sim"


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
#' @format ## `zhang_sig_activity`
#' A dataframe with 39 rows and 213 columns:
#' \describe{
#' \item{Subject}{Unique ID for each subject}
#' ...
#' \item{passive_smoking}{Subject smoking status ("Non-smoker" or "Passive smoker")}
#' ...
#' \item{SBS1, ..., SBS40}{Relative abundance of mutational signatures in each sample. These columns sum to 1 for each row.}
#' }
#' @source <https://doi.org/10.1038/s41588-021-00920-0>
"zhang_sig_activity"


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

#' Single base substitution signature attributions found in esophageal squamous cell carcinoma across countries with varying incidence.
#'
#' Data from: Moody et al. Nat Genet (2021)
"ESCC_countries_SBS"

#' Reference COSMIC Single base substitution signatures found in esophageal squamous cell carcinoma across countries with varying incidence.
#'
#' Data from: Moody et al. Nat Genet (2021)
"ESCC_countries_SBS.refs"

"Sherlock_LCINS.metadata"

#' Single base substitution signature attributions found in Lung cancer in never smokers.
#'
#' Data from: Zhang et al. Nat Genet (2021)
"Sherlock_LCINS_SBS"

#' Reference COSMIC Single base substitution signatures found in Lung cancer in never smokers.
#'
#' Data from: Zhang et al. Nat Genet (2021)
"Sherlock_LCINS_SBS.refs"

#' Single base substitution signature attributions found in mice exposed to 20 known or suspected carcinogens.
#'
#' Data from: Riva et al. Nat Genet (2020)
"carcinogens_mice_SBS"


#' Reference Single base substitution signatures found in mice exposed to 20 known or suspected carcinogens.
#'
#' Data from: Riva et al. Nat Genet (2020)
"carcinogens_mice_SBS.refs"


#' Single base substitution signature attributions found in papillary thyroid carcinoma from Chernobyl incident survivors.
#'
#' Data from: Morton et al. Science (2021)
"radiation_sigs_morton"

#' Cosine similarity matrix of Single base substitution signatures found in papillary thyroid carcinoma from Chernobyl incident survivors.
#'
#' Data from: Morton et al. Science (2021)
"radiation_sigs_morton_cossim"

#' Cosine similarity matrix of Single base substitution signatures found in lung adenocarcinomas in Asia.
#'
#' Data from: Chen et al. Nat Genet (2020)
"smoker_sigs_chen_cossim"


#' Supplementary table from Moody et al. Nat Genet 2021.
#'
#' Data from: Moody et al. Nat Genet (2021)
"tab15"


#' Driver mutations from the intogen website for gene EGFR in lung adenocarcinomas
#'
#' Data from: Intogen database
"Intogen_EGFR_LUAD"

#' Driver mutations from the intogen website for gene TP53 in lung adenocarcinomas
#'
#' Data from: Intogen database
"TP53_drivers_intogen_LUAD"
