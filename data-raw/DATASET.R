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

  "mSBS42" = "#FD6CFD", #Ultra pink

  mSBS_N1 = "#4D4D4D"
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

# Get Sherlock data
Sherlock_LCINS_SBS = read_tsv("/home/alcalan/evolution/data/signatures/non-smoker_lung_cancer_genomics/sherlock_232_SBS_Signature_Ratio.txt")
Sherlock.muts     = read_xlsx("/home/alcalan/evolution/data/signatures/non-smoker_lung_cancer_genomics/41588_2021_920_MOESM4_ESM.xlsx",sheet=4,skip=1)
Sherlock_LCINS_SBS.refs = read_xlsx("/home/alcalan/evolution/data/signatures/non-smoker_lung_cancer_genomics/41588_2021_920_MOESM4_ESM.xlsx",sheet=7,skip=2) %>% arrange(Type,Subtype)
Sherlock_LCINS.metadata = read_xlsx("/home/alcalan/evolution/data/signatures/non-smoker_lung_cancer_genomics/41588_2021_920_MOESM4_ESM.xlsx",sheet=1,skip=1)

Sherlock_LCINS_SBS.refs = Sherlock_LCINS_SBS.refs[,colnames(Sherlock_LCINS_SBS)[-1]]
all(colnames(Sherlock_LCINS_SBS.refs)==colnames(Sherlock_LCINS_SBS)[-(1)])

Sherlock_LCINS_SBS.refs$Type = MutType$Type
Sherlock_LCINS_SBS.refs$Subtype = MutType$Subtype

usethis::use_data(Sherlock_LCINS_SBS,overwrite = TRUE)
usethis::use_data(Sherlock_LCINS_SBS.refs,overwrite = TRUE)
usethis::use_data(Sherlock_LCINS.metadata,overwrite = TRUE)

# code to prepare mice carcinogen exposure dataset
mice_metadata = readxl::read_xlsx("../data/mice_carcinogens/41588_2020_692_MOESM2_ESM.xlsx",sheet = 2,skip=3)
mice_metadata = mice_metadata[-1,] %>% rename(Sample=`Sample name in the manuscript`)
mice_metadata$Sample = str_replace(mice_metadata$Sample,"VANDIUM","VANADIUM") # correct name

mexposure = readRDS("../data/mouse-mutatation-signatures/figure1/mexposure.rds")
mexposure.tib = bind_cols(Sample=str_replace(stringr::str_replace_all( rownames(mexposure) ," ","_"),"STOMACH","FORESTOMACH"), as_tibble(mexposure))

# merge
mexposure.tib$Sample[!mexposure.tib$Sample %in% mice_metadata$Sample] # all here

carcinogens_mice_SBS = left_join(mexposure.tib,mice_metadata)

# clean dose
carcinogens_mice_SBS$dose_numeric = as.numeric( str_remove(mutsig_carcinogens_mice_SBS$dose," ppm| mg/m3| MG/KG| PPM| MG/L| mg/L| mg/kg| m.9ful|mg/m3| mg/l| MG/M3") )

# reorder columns
carcinogens_mice_SBS = mutsig_carcinogens_mice_SBS %>% dplyr::relocate(.after = Sample, mSBS1,mSBS5,mSBS12,mSBS17,mSBS18,mSBS19,mSBS40, mSBS42, mSBS_N1, mSBS_N2, mSBS_N3)

# save
usethis::use_data(carcinogens_mice_SBS,overwrite = TRUE)

## get references for mSBS
carcinogens_mice_SBS.refs = readRDS("../data/mouse-mutatation-signatures/figure1/mSBSs.rds")

barplot( mutsig_carcinogens_mice_SBS.refs[,1] ) # to check order. Seems to be the same as in fig 1d (C>A, C>G, C>T, T>A, etc)
#carcinogens_mice_SBS.refs = bind_cols( Type=Sherlock_LCINS_SBS.refs$Type,  Subtype=Sherlock.refs$Subtype, as_tibble(carcinogens_mice_SBS.refs))

carcinogens_mice_SBS.refs = mutsig_carcinogens_mice_SBS.refs[,sigs_mice.S]

carcinogens_mice_SBS.refs$Type = MutType$Type
carcinogens_mice_SBS.refs$Subtype = MutType$Subtype

## save
usethis::use_data(carcinogens_mice_SBS.refs,overwrite = TRUE)

## read mutational spectra data
mutprof_carcinogens_mice_SBS = read_tsv("../data/signatures/COSMIC/SPinput_split_mice_carcinogens/output/SBS/Mice_Carcinogens.SBS96.all") %>%
  mutate(Type =  str_extract(MutationType,"[ACGT]>[ACGT]"),
         Subtype = paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7)) ) %>%
  arrange(Type,Subtype) %>% relocate(Type,Subtype,.after= MutationType)

### save
usethis::use_data(mutprof_carcinogens_mice_SBS,overwrite = TRUE)

## transcripts
### get sequences
Hras.trans = read_tsv("../data/Ensembl/Mus_musculus_Hras_ENSMUST00000026572_11_sequence.fa")
Hras.trans = unlist( str_split(Hras.trans$`>Hras-201 ENSMUSE00000480008 exon:protein_coding`[!str_detect(Hras.trans$`>Hras-201 ENSMUSE00000480008 exon:protein_coding`,">")],"") )

Kras.trans = read_tsv("../data/Ensembl/Mus_musculus_Kras_ENSMUST00000111710_8_sequence.fa")
Kras.trans = unlist( str_split(Kras.trans$`>Kras-202 ENSMUSE00000380468 exon:protein_coding`[!str_detect(Kras.trans$`>Kras-202 ENSMUSE00000380468 exon:protein_coding`,">")],"") )

Fgfr2.trans = read_tsv("../data/Ensembl/Mus_musculus_Fgfr2_ENSMUST00000122054_8_sequence.fa")
Fgfr2.trans = unlist( str_split(Fgfr2.trans$`>Fgfr2-215 ENSMUSE00000721980 exon:protein_coding`[!str_detect(Fgfr2.trans$`>Fgfr2-215 ENSMUSE00000721980 exon:protein_coding`,">")],"") )

Braf.trans = read_tsv("../data/Ensembl/Mus_musculus_Braf_ENSMUST00000002487_15_sequence.fa")
Braf.trans = unlist( str_split(Braf.trans$`>Braf-201 ENSMUSE00000889921 exon:protein_coding`[!str_detect(Braf.trans$`>Braf-201 ENSMUSE00000889921 exon:protein_coding`,">")],"") )

compute_SBS96_context <- function(trans){
### create contexts
trans.contexts  = sapply(2:(length(trans)-1),  function(i) paste0(trans[i+(-1):1],collapse = "") )

### Convert to SBS96
reverse <- function(X){
  res = str_split(X,"")[[1]]
  res = sapply(res, function(x){
    if(x=="A") return("T")
    if(x=="T") return("A")
    if(x=="C") return("G")
    if(x=="G") return("C")
  })
  return(paste0(res,collapse = ""))
}

### apply when necessary
trans.contexts.96 = sapply(2:(length(trans)-1), function(i){
  if(trans[i] %in% c("G","A") ){
    reverse(trans.contexts[i-1])
  }else{
    trans.contexts[i-1]
    }
  } )

return( table(trans.contexts.96)[MutType$Subtype] )
}

Braf.trans.SBS96.context  = compute_SBS96_context(Braf.trans)
Kras.trans.SBS96.context  = compute_SBS96_context(Kras.trans)
Hras.trans.SBS96.context  = compute_SBS96_context(Hras.trans)
Fgfr2.trans.SBS96.context = compute_SBS96_context(Fgfr2.trans)

# code to create mice mutsig files
#Mice_Carcinogens.SBS96.input0 = read_tsv("../data/signatures/mice_carcinogens/41588_2020_692_MOESM2_ESM.xlsx")
Mice_Carcinogens.SBS96.input = read_tsv("../data/signatures/COSMIC/SPinput_split_mice_carcinogens/Mice_carcinogens_driver_SPformat.txt")
Mice_Carcinogens.SBS96 = read_tsv("../data/signatures/COSMIC/SPinput_split_mice_carcinogens/output/SBS/Mice_Carcinogens.SBS96.all") # contexts looks ok compared to UCSC

# drivers
## reduced list
drivers = read_xlsx("../data/mice_carcinogens/41588_2020_692_MOESM2_ESM.xlsx",sheet = 12,skip=2)
# whole table of Fig. 4
drivers_fig4  = read_tsv("../data/mouse-mutatation-signatures/figure4/snvsindrivers.txt",col_names = NA) # subset only to the ones above 3% as in Fig4a
drivers_fig4a = read_tsv("../data/mouse-mutatation-signatures/figure4/drivers_lung.txt",col_names = NA)
drivers_fig4b = read_tsv("../data/mouse-mutatation-signatures/figure4/drivers_liver.txt",col_names = NA)

drivers_fig4.lung = drivers_fig4a %>% filter(X6%in%c("Kras","Fgfr2","Braf","Ctnnb1","Zfhx3","Trp53","Epha3","Cnot3"))
drivers_fig4.liver = drivers_fig4b %>% filter(X6%in%c("Hras","Ctnnb1","Egfr","Lrp1b","Braf","Kras"))

drivers_fig4.lung = drivers_fig4.lung %>% mutate(mutation_ID=paste(str_replace_all(X1," ","_"),X2,X3,X4,X5,sep="_"))
drivers_fig4.liver = drivers_fig4.liver %>% mutate(mutation_ID=paste(str_replace_all(X1," ","_"),X2,X3,X4,X5,sep="_"))

# read sig outputs
mut_profs_all_drivers = read_tsv("../data/signatures/COSMIC/SPinput_split_mice_carcinogens/output/SBS/Mice_Carcinogens.SBS96.all") %>% 
  mutate(Type=str_extract(MutationType,"[ATCG]>[ATCG]"),
         Subtype=paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7)) ) %>% dplyr::arrange(Type,Subtype)


mut_profs_all_drivers.t = t(mut_profs_all_drivers[,-c(1,ncol(mut_profs_all_drivers)-1,ncol(mut_profs_all_drivers))])
colnames(mut_profs_all_drivers.t) = mut_profs_all_drivers$MutationType
mut_profs_all_drivers.t = bind_cols(mutation_ID= rownames(mut_profs_all_drivers.t), as_tibble(mut_profs_all_drivers.t) )

drivers_fig4.lung.spectra  = left_join(drivers_fig4.lung,mut_profs_all_drivers.t)
drivers_fig4.liver.spectra = left_join(drivers_fig4.liver,mut_profs_all_drivers.t)

mice_lung_drivers_SBS = tibble(Fgfr2 = colSums(drivers_fig4.lung.spectra[drivers_fig4.lung.spectra$X6=="Fgfr2",1:96+10]) ,
                               Kras = colSums(drivers_fig4.lung.spectra[drivers_fig4.lung.spectra$X6=="Kras",1:96+10])
                               )

mice_liver_drivers_SBS = tibble(Hras = colSums(drivers_fig4.liver.spectra[drivers_fig4.liver.spectra$X6=="Hras",1:96+10]) ,
                               Braf = colSums(drivers_fig4.liver.spectra[drivers_fig4.liver.spectra$X6=="Braf",1:96+10])
)

mice_lung_drivers_SBS$Fgfr2 = mice_lung_drivers_SBS$Fgfr2/Fgfr2.trans.SBS96.context
mice_lung_drivers_SBS$Kras = mice_lung_drivers_SBS$Fgfr2/Kras.trans.SBS96.context
mice_liver_drivers_SBS$Hras = mice_liver_drivers_SBS$Hras/Hras.trans.SBS96.context
mice_liver_drivers_SBS$Braf = mice_liver_drivers_SBS$Braf/Braf.trans.SBS96.context

usethis::use_data(mice_liver_drivers_SBS,overwrite = T)
usethis::use_data(mice_lung_drivers_SBS,overwrite = T)
#

sigs_mice.S = colnames(mutsig_carcinogens_mice_SBS[,2:12])

Riva_SBS       = as.matrix(mutsig_carcinogens_mice_SBS[,2:12])
rownames(Riva_SBS) = mutsig_carcinogens_mice_SBS$Sample
Riva.refs = as.matrix(mutsig_carcinogens_mice_SBS.refs[,sigs_mice.S])
colnames(Riva.refs) = sigs_mice.S
rownames(Riva.refs) = paste0(mutsig_carcinogens_mice_SBS.refs$Type,",",mutsig_carcinogens_mice_SBS.refs$Subtype)

Riva_MutProf = Riva_SBS%*%t(Riva.refs)

drivers = drivers %>% mutate(Tissue = str_extract(sample,"[A-Z]+"), 
                             Chemical = str_remove_all(str_extract(sample,"_[A-Z_]+"),"^_|_$"), 
                             Type=str_extract(context,"[ATCG]>[ATCG]"),
                             Subtype=paste0(str_sub(context,1,1),str_sub(context,3,3),str_sub(context,7,7)) )

drivers.spectrum.lung.mice = drivers[!duplicated(drivers[,2:5]),] %>% filter(Tissue=="LUNG") %>% group_by(gene,Type,Subtype) %>% summarize(n=n())
drivers.spectrum.lung.mice.l = lapply(unique(drivers.spectrum.lung.mice$gene), function(x)  left_join(MutType,drivers.spectrum.lung.mice %>% filter(gene==x)) )
for(i in 1:length(drivers.spectrum.lung.mice.l)) drivers.spectrum.lung.mice.l[[i]]$n[is.na(drivers.spectrum.lung.mice.l[[i]]$n)] = 0
names(drivers.spectrum.lung.mice.l) = unique(drivers.spectrum.lung.mice$gene)

drivers.spectrum.liver.mice = drivers[!duplicated(drivers[,2:5]),] %>% filter(Tissue=="LIVER") %>% group_by(gene,Type,Subtype) %>% summarize(n=n())
drivers.spectrum.liver.mice.l = lapply(unique(drivers.spectrum.liver.mice$gene), function(x)  left_join(MutType,drivers.spectrum.liver.mice %>% filter(gene==x)) )
for(i in 1:length(drivers.spectrum.liver.mice.l)) drivers.spectrum.liver.mice.l[[i]]$n[is.na(drivers.spectrum.liver.mice.l[[i]]$n)] = 0
names(drivers.spectrum.liver.mice.l) = unique(drivers.spectrum.liver.mice$gene)

drivers.spectrum.lung.mice.l$Kras$n/Kras.trans.SBS96.context

mice_lung_drivers_SBS
mice_liver_drivers_SBS

# ESCC signatures
load("data/tab15.rda")

ESCC_countries_SBS = tab15[,!str_detect(colnames(tab15),"ID|DBS")]
# renormalize
ESCC_countries_SBS[,-c(1:3)] = sweep(ESCC_countries_SBS[,-c(1:3)] , 1, rowSums(ESCC_countries_SBS[,-c(1:3)]), "/")

all_sigs = colnames(ESCC_countries_SBS)[-(1:3)]

escc_stats.SBS = sigvar::sigvar(sig_activity = ESCC_countries_SBS, K = length(all_sigs), group = "Country"#, S = ESCC_sim
) %>%
  mutate(incidence = ifelse(Country %in% c("UK", "Japan", "Brazil"),
                            "Low incidence", "High incidence"))

# References
SBS_ref = read_tsv("../data/signatures/COSMIC/COSMIC_v3.2_SBS_GRCh37.txt")
SBS_denovo = readxl::read_xlsx("../data/signatures/Mutographs_signatures/MoodyEtAl2021_suptables.xlsx",
                               sheet = "Supplementary Table 2", skip = 2)
SBS_all = left_join(SBS_ref, SBS_denovo %>%
                      mutate(MutationsType = stringr::str_remove(MutationsType, "N\\:")),
                    by = c("Type" = "MutationsType")) %>%
  mutate(MutationType=Type,
         Type = str_extract(MutationType,"[ACGT]>[ACGT]"),
         Subtype= paste0(str_sub(MutationType,1,1),
                         str_sub(MutationType,3,3),
                         str_sub(MutationType,7,7)))

ESCC_countries_SBS.refs = left_join(MutType,SBS_all)

sigs_ESCC.S = colnames(ESCC_countries_SBS)[-(1:3)]
ESCC_countries_SBS.refs = ESCC_countries_SBS.refs[,sigs_ESCC.S]

ESCC_countries_SBS.refs$Type = MutType$Type
ESCC_countries_SBS.refs$Subtype = MutType$Subtype

# save data
usethis::use_data(ESCC_countries_SBS,overwrite = T)
usethis::use_data(ESCC_countries_SBS.refs,overwrite = T)

# ESCC drivers
ESCC_drivers = read_xlsx("../data/signatures/Mutographs_signatures/MoodyEtAl2021_suptables.xlsx",sheet=12,skip=2)
ESCC_drivers.SBS = ESCC_drivers %>% filter(Type=="Sub")

sort(table(ESCC_drivers.SBS$Gene),decreasing=T) # TP53 overwhelmingly greater

# write in SP format: Project	Sample	ID	Genome	mut_type	chrom	pos_start	pos_end	ref	alt	Type
ESCC_drivers.SBS.SPformat = ESCC_drivers.SBS %>%
  mutate(Project="ESCC",Sample= paste(Sample,Chrom,Pos,Ref,Alt,sep="_"), ID=Gene,Genome="GRCh37",mut_type=Type,chrom=Chrom,pos_start=Pos,	pos_end=Pos,
         ref=Ref,	alt=Alt,	Type="SOMATIC") %>%
  dplyr::select(Project	,Sample,	ID,	Genome,	mut_type,	chrom,	pos_start,	pos_end,	ref,	alt,Type)

write_tsv(ESCC_drivers.SBS.SPformat,file="../data/signatures/COSMIC/SPinput_ESCC/ESCC_drivers_SPformat.txt")

# read driver spectra
ESCC_drivers.SBS.spectra = read_tsv("../data/signatures/COSMIC/SPinput_ESCC/output/SBS/COSMIC_EGFR.SBS96.all")

ESCC_drivers.SBS.spectra = left_join(MutType , ESCC_drivers.SBS.spectra %>% mutate(Type = str_extract(MutationType,"[ACGT]>[ACGT]"),
                                                        Subtype= paste0(str_sub(MutationType,1,1),
                                                                        str_sub(MutationType,3,3),
                                                                        str_sub(MutationType,7,7))
                                                        ) )
ESCC_drivers.SBS.spectra[,-(1:3)] = ESCC_drivers.SBS.spectra[,-(1:3)][sort(colnames(ESCC_drivers.SBS.spectra)[-(1:3)])]
ESCC_drivers.SBS.SPformat = ESCC_drivers.SBS.SPformat %>% arrange(Sample)
#check match
all(ESCC_drivers.SBS.SPformat$Sample==colnames(ESCC_drivers.SBS.spectra)[-(1:3)])
## get TP53 spectrum
ESCC_TP53driverpos.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="TP53" &
                                                                                                                         !duplicated(ESCC_drivers.SBS.SPformat[,6:10]))+3]))

#ESCC_TP53drivers.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="TP53")+3]))
#ESCC_CDKN2Adrivers.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="CDKN2A")+3]))
ESCC_CDKN2Adriverpos.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="CDKN2A" &
                                                                                                                         !duplicated(ESCC_drivers.SBS.SPformat[,6:10]))+3]))

#ESCC_PIK3CAdrivers.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="PIK3CA")+3]))
#ESCC_PIK3CAdriverpos.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="PIK3CA" &
#                                                                                                                           !duplicated(ESCC_drivers.SBS.SPformat[,6:10]))+3]))

ESCC_KMT2Ddrivers.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="KMT2D")+3]))
ESCC_KMT2Ddriverpos.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="KMT2D" &
                                                                                                                           !duplicated(ESCC_drivers.SBS.SPformat[,6:10]))+3]))


ESCC_EP300drivers.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="EP300")+3]))
ESCC_EP300driverpos.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="EP300" &
                                                                                                                           !duplicated(ESCC_drivers.SBS.SPformat[,6:10]))+3]))

ESCC_NFE2L2drivers.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="NFE2L2")+3]))
ESCC_NFE2L2driverpos.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="NFE2L2" &
                                                                                                                           !duplicated(ESCC_drivers.SBS.SPformat[,6:10]))+3]))

ESCC_NOTCH1drivers.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="NOTCH1")+3]))
ESCC_NOTCH1driverpos.SBS.spectra = bind_cols(ESCC_drivers.SBS.spectra[,c(1:3)],n=rowSums(ESCC_drivers.SBS.spectra[,which(ESCC_drivers.SBS.SPformat$ID=="NOTCH1" &
                                                                                                                           !duplicated(ESCC_drivers.SBS.SPformat[,6:10]))+3]))

# Compute contexts
TP53.trans = unlist(str_split(read_tsv("../data/signatures/COSMIC/Homo_sapiens_TP53_ENST00000269305_9_sequence.fa",comment = ">",col_names = NA)[[1]],""))
TP53.trans.SBS96.context  = compute_SBS96_context(TP53.trans)

EGFR.trans = unlist(str_split(read_tsv("../data/signatures/COSMIC/Homo_sapiens_EGFR_ENST00000275493_7_sequence.fa",comment = ">",col_names = NA)[[1]],""))
EGFR.trans.SBS96.context  = compute_SBS96_context(EGFR.trans)

CDKN2A.trans = unlist(str_split(read_tsv("../data/signatures/COSMIC/Homo_sapiens_CDKN2A_ENST00000304494_10_sequence.fa",comment = ">",col_names = NA)[[1]],""))
CDKN2A.trans.SBS96.context  = compute_SBS96_context(CDKN2A.trans)

PIK3CA.trans = unlist(str_split(read_tsv("../data/signatures/COSMIC/Homo_sapiens_PIK3CA_ENST00000263967_4_sequence.fa",comment = ">",col_names = NA)[[1]],""))
PIK3CA.trans.SBS96.context  = compute_SBS96_context(PIK3CA.trans)

KMT2D.trans = unlist(str_split(read_tsv("../data/signatures/COSMIC/Homo_sapiens_KMT2D_ENST00000301067_12_sequence.fa",comment = ">",col_names = NA)[[1]],""))
KMT2D.trans.SBS96.context  = compute_SBS96_context(KMT2D.trans)

EP300.trans = unlist(str_split(read_tsv("../data/signatures/COSMIC/Homo_sapiens_EP300_ENST00000263253_9_sequence.fa",comment = ">",col_names = NA)[[1]],""))
EP300.trans.SBS96.context  = compute_SBS96_context(EP300.trans)

NFE2L2.trans = unlist(str_split(read_tsv("../data/signatures/COSMIC/Homo_sapiens_NFE2L2_ENST00000397062_8_sequence.fa",comment = ">",col_names = NA)[[1]],""))
NFE2L2.trans.SBS96.context  = compute_SBS96_context(NFE2L2.trans)

NOTCH1.trans = unlist(str_split(read_tsv("../data/signatures/COSMIC/Homo_sapiens_NOTCH1_ENST00000651671_1_sequence.fa",comment = ">",col_names = NA)[[1]],""))
NOTCH1.trans.SBS96.context  = compute_SBS96_context(NOTCH1.trans)

# driver spectra
ESCC_drivers_SBS = tibble(Type=ESCC_TP53driverpos.SBS.spectra$Type,
                          Subtype=ESCC_TP53driverpos.SBS.spectra$Subtype,
  TP53=ESCC_TP53driverpos.SBS.spectra$n/TP53.trans.SBS96.context,
  CDKN2A=ESCC_CDKN2Adriverpos.SBS.spectra$n/TP53.trans.SBS96.context,
  EP300=ESCC_EP300driverpos.SBS.spectra$n/TP53.trans.SBS96.context,
  KMT2D=ESCC_KMT2Ddriverpos.SBS.spectra$n/TP53.trans.SBS96.context,
  NFE2L2=ESCC_NFE2L2driverpos.SBS.spectra$n/TP53.trans.SBS96.context,
  NOTCH1=ESCC_NOTCH1driverpos.SBS.spectra$n/TP53.trans.SBS96.context,
  PIK3CA=ESCC_PIK3CAdriverpos.SBS.spectra$n/TP53.trans.SBS96.context )

usethis::use_data(ESCC_drivers_SBS,overwrite = T)

## LUAD drivers intogen
EGFR_drivers_intogen_SBS96 = read_tsv("../data/signatures/COSMIC/SPinput_split/output/SBS/COSMIC_EGFR.SBS96.all")
EGFR_drivers_intogen_SBS96 = tibble( MutationType=EGFR_drivers_intogen_SBS96$MutationType, EGFR=rowSums(EGFR_drivers_intogen_SBS96[,-1]) ) %>% 
  mutate(Type=str_extract(MutationType,"[ATCG]>[ATCG]"),
         Subtype=paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7)) ) %>% dplyr::arrange(Type,Subtype)
TP53_drivers_intogen_SBS96 = read_tsv("../data/signatures/COSMIC/SPinput_TP53_split/output/SBS/COSMIC_EGFR.SBS96.all")
TP53_drivers_intogen_SBS96 = tibble( MutationType=TP53_drivers_intogen_SBS96$MutationType, TP53=rowSums(TP53_drivers_intogen_SBS96[,-1]) )%>% 
  mutate(Type=str_extract(MutationType,"[ATCG]>[ATCG]"),
         Subtype=paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7)) ) %>% dplyr::arrange(Type,Subtype)

LUAD_drivers_SBS = left_join(left_join(MutType,EGFR_drivers_intogen_SBS96),TP53_drivers_intogen_SBS96)

LUAD_drivers_SBS$EGFR = LUAD_drivers_SBS$EGFR/EGFR.trans.SBS96.context
LUAD_drivers_SBS$TP53 = LUAD_drivers_SBS$TP53/TP53.trans.SBS96.context

LUAD_drivers_SBS = LUAD_drivers_SBS[,-3]
usethis::use_data(LUAD_drivers_SBS,overwrite = T)
