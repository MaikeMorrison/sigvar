# code with signature colors
sbs_palette = c(
  # shades of royal blue
  "SBS40" = "#2874A6",
  "SBS40a" = "#81A1C8",
  "SBS40b" = "#0E4384",
  "SBS40c" = "#2874A6",

  "SBS5" = "#FFC642", # mustard
  "SBS1" = "#C0392B", # brick
  "SBS13" = "#388E3C", # lizard. contrasts with 18
  "SBS18" = "#FFA726", # marigold -
  "SBS2" = "#223859", # navy
  "SBS4" = "#8BC34A", # green - seafoam81C784 contrasts with 40a/b, 5
  "SBS3" = "#F8C471", # soft pumpkin
  "SBS12" = "#633974", #amethyst

  "SBS22" = "#D35400",  # dark terra cotta
  "SBS22b" = "#DC7633", # light terra cotta
  "SBS22a" = "#D35400", # dark terra cotta. contrasts with 40a/b

  "SBS17a" = "#B2EBF2",  # light turquoise
  "SBS17b" = "#00ACC1", # dark turquoise

  "SBS19" = "#00A86B",  # Jade green

  "SBS7a" = "#EC407A", # magenta
  "SBS16" = "#8E44AD", # purple
  "SBS8" = "#F9E79F", # soft yellow
  "SBS9" = "#ABB2B9", # grey
  "SBS29" = "#F5B7B1", # soft pink
  "SBS41" = "#B39DDB", #periwinkle

  "SBS42" = "#FD6CFD", #Ultra pink

  ### mice signatures
  "mSBS40" = "#2874A6",
  "mSBS5" = "#FFC642", # mustard
  "mSBS1" = "#C0392B", # brick
  "mSBS13" = "#388E3C", # lizard. contrasts with 18
  "mSBS18" = "#FFA726", # marigold -
  "mSBS2" = "#223859", # navy
  "mSBS4" = "#8BC34A", # green - seafoam81C784 contrasts with 40a/b, 5
  "mSBS3" = "#F8C471", # soft pumpkin
  "mSBS12" = "#633974", #amethyst

  "mSBS22" = "#D35400",  # dark terra cotta
  "mSBS22b" = "#DC7633", # light terra cotta
  "mSBS22a" = "#D35400", # dark terra cotta. contrasts with 40a/b

  "mSBS17" = "#B2EBF2",  # light turquoise

  "mSBS19" = "#00A86B",  # Jade green

  "mSBS7a" = "#EC407A", # magenta
  "mSBS16" = "#8E44AD", # purple
  "mSBS8" = "#F9E79F", # soft yellow
  "mSBS9" = "#ABB2B9", # grey
  "mSBS29" = "#F5B7B1", # soft pink
  "mSBS41" = "#B39DDB", #periwinkle

  "mSBS42" = "#FD6CFD" #Ultra pink
)

sbs_palette = sbs_palette[order(names(sbs_palette))]
usethis::use_data(sbs_palette,overwrite = T)

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
