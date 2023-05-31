# code to prepare PCAWG dataset
PCAWG_SigProfiler_COSMIC_SBS = read_csv("../data/signatures/PCAWG_signatures/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
PCAWG_SigProfiler_COSMIC_DBS = read_csv("../data/signatures/PCAWG_signatures/PCAWG_sigProfiler_DBS_signatures_in_samples.csv")
PCAWG_SigProfiler_COSMIC_ID  = read_csv("../data/signatures/PCAWG_signatures/PCAWG_SigProfiler_ID_signatures_in_samples.csv")

usethis::use_data(PCAWG_SigProfiler_COSMIC_SBS)
usethis::use_data(PCAWG_SigProfiler_COSMIC_DBS)
usethis::use_data(PCAWG_SigProfiler_COSMIC_ID)

SPfolder = system.file("extdata", "example_SigProfiler_results", package = "sigFAVA")
internal_for_import_test = import_SigProfiler(SPfolder)
usethis::use_data(internal_for_import_test, internal = TRUE)

# code to prepare mice carcinogen exposure dataset
mice_metadata = readxl::read_xlsx("../data/mice_carcinogens/41588_2020_692_MOESM2_ESM.xlsx",sheet = 2,skip=3)
mice_metadata = mice_metadata[-1,] %>% rename(Sample=`Sample name in the manuscript`)
mice_metadata$Sample = str_replace(mice_metadata$Sample,"VANDIUM","VANADIUM") # correct name

mexposure = readRDS("../data/mouse-mutatation-signatures/figure1/mexposure.rds")
mexposure.tib = bind_cols(Sample=str_replace(stringr::str_replace_all( rownames(mexposure) ," ","_"),"STOMACH","FORESTOMACH"), as_tibble(mexposure))

# merge
mexposure.tib$Sample[!mexposure.tib$Sample %in% mice_metadata$Sample] # all here

mutsig_carcinogens_mice_SBS = left_join(mexposure.tib,mice_metadata)

# clean dose
mutsig_carcinogens_mice_SBS$dose_numeric = as.numeric( str_remove(mutsig_carcinogens_mice$dose," ppm| mg/m3| MG/KG| PPM| MG/L| mg/L| mg/kg| m.9ful|mg/m3| mg/l| MG/M3") )

# reorder columns
mutsig_carcinogens_mice_SBS = mutsig_carcinogens_mice %>% dplyr::relocate(.after = Sample, mSBS1,mSBS5,mSBS12,mSBS17,mSBS18,mSBS19,mSBS40, mSBS42, mSBS_N1, mSBS_N2, mSBS_N3)

# save
usethis::use_data(mutsig_carcinogens_mice_SBS,overwrite = TRUE)
