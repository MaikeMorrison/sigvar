## code to prepare dataset
PCAWG_SigProfiler_COSMIC_SBS = read_csv("../data/signatures/PCAWG_signatures/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
PCAWG_SigProfiler_COSMIC_DBS = read_csv("../data/signatures/PCAWG_signatures/PCAWG_sigProfiler_DBS_signatures_in_samples.csv")
PCAWG_SigProfiler_COSMIC_ID  = read_csv("../data/signatures/PCAWG_signatures/PCAWG_SigProfiler_ID_signatures_in_samples.csv")

usethis::use_data(PCAWG_SigProfiler_COSMIC_SBS)
usethis::use_data(PCAWG_SigProfiler_COSMIC_DBS)
usethis::use_data(PCAWG_SigProfiler_COSMIC_ID)

SPfolder = system.file("extdata", "example_SigProfiler_results", package = "sigFAVA")
internal_for_import_test = import_SigProfiler(SPfolder)
usethis::use_data(internal_for_import_test, internal = TRUE)
