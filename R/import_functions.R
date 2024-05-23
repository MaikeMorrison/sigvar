#' import_SigProfiler
#'
#' @param folder a folder containing SigProfilerExtractor results
#'
#' @return A list of tibbles containing the SigProfiler activity matrices
#' @export
#'
#' @examples
#' #Import dummy SigProfilerExtractor output results
#' SPfolder = system.file("extdata", "SP", package = "sigvar")
#' Qlist = import_SigProfiler(SPfolder)
#' print(Qlist)
#' for(i in 1:length(Qlist)) Qlist[[i]][,-1] = sweep(Qlist[[i]][,-1],1,rowSums(Qlist[[i]][,-1]),"/")
#'
#' #Plot imported signatures
#' plot_dots(Qlist[[1]])
#' @importFrom readr read_tsv
import_SigProfiler <- function(folder="."){
  # find activity files within the folder
  input_files        = list.files(folder,pattern = "Activities[_refit]*.txt",recursive = T,full.names = T)
  # get name (de novo or COSMIC)
  input_files.names  = sapply( list.files(folder,pattern = "Activities[_refit]*.txt",recursive = T,full.names = F), function(x){strsplit(x ,"/")[[1]][3]} )
  # keep only suggested solutions instead of all NMF solutions
  input_files.tokeep = grep(input_files, pattern = "Suggested" ,value = F)
  # read files
  Qlist = lapply( input_files[input_files.tokeep], read_tsv)
  # rename list entries
  names(Qlist) = input_files.names[input_files.tokeep]
  return(Qlist)
}


#' get_SBS96_spectrum
#'
#' @param transcript the ensembl ID of the transcript
#' @param organism the name of the organism associated with the transcript (species)
#'
#' @return A vector of size 96 containing the SBS spectrum of the transcript
#' @export
#'
#' @examples
#' #Run on
#' spectrum = get_SBS96_spectrum(transcript = "ENST00000269305.9")
#' print(spectrum)
#' @importFrom GenomicFeatures makeTxDbFromBiomart
#' @importFrom GenomicFeatures exons
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom Biostrings getSeq
#' @importFrom rtracklayer chrom
#' @importFrom rtracklayer start
#' @importFrom rtracklayer end
#' @importFrom rtracklayer strand
#' @importFrom Biostrings trinucleotideFrequency
#' @importFrom Biostrings complement
#' @importFrom Biostrings reverse
#' @importFrom Biostrings DNAStringSet
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
get_SBS96_spectrum = function(transcript="ENST00000269305.9",organism = "Homo sapiens"){
  if(sub("_| ","",tolower(organism))=="homosapiens") txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  if(!exists("txdb")) stop("Organism not found. Valid answers are homo sapiens")
  SBS32_Subtypes = c("ACA", "ACC", "ACG", "ACT",
                     "CCA", "CCC", "CCG", "CCT",
                     "GCA", "GCC", "GCG", "GCT",
                     "TCA", "TCC", "TCG", "TCT",
                     "ATA", "ATC", "ATG", "ATT",
                     "CTA", "CTC", "CTG", "CTT",
                     "GTA", "GTC", "GTG", "GTT",
                     "TTA", "TTC", "TTG", "TTT")
  SBS96_Subtypes = c("ACA", "ACC", "ACG", "ACT","CCA", "CCC", "CCG", "CCT","GCA", "GCC", "GCG", "GCT","TCA", "TCC", "TCG", "TCT",
                     "ACA", "ACC", "ACG", "ACT","CCA", "CCC", "CCG", "CCT","GCA", "GCC", "GCG", "GCT","TCA", "TCC", "TCG", "TCT",
                     "ACA", "ACC", "ACG", "ACT","CCA", "CCC", "CCG", "CCT","GCA", "GCC", "GCG", "GCT","TCA", "TCC", "TCG", "TCT",
                     "ATA", "ATC", "ATG", "ATT","CTA", "CTC", "CTG", "CTT","GTA", "GTC", "GTG", "GTT","TTA", "TTC", "TTG", "TTT",
                     "ATA", "ATC", "ATG", "ATT","CTA", "CTC", "CTG", "CTT","GTA", "GTC", "GTG", "GTT","TTA", "TTC", "TTG", "TTT",
                     "ATA", "ATC", "ATG", "ATT","CTA", "CTC", "CTG", "CTT","GTA", "GTC", "GTG", "GTT","TTA", "TTC", "TTG", "TTT")
  # retrieve transcript sequence + neighboring nucleotides from DB
  hg38_tr = GenomicFeatures::exons(txdb,filter=list(tx_name=transcript))
  GenomeInfoDb::seqlevelsStyle(hg38_tr) <- GenomeInfoDb::seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  hg38_seq <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, names=rtracklayer::chrom(hg38_tr) ,
                                 start=rtracklayer::start(hg38_tr)-1,end=rtracklayer::end(hg38_tr)+1,
                                 strand=rtracklayer::strand(hg38_tr))
  # compute spectrum
  spectrum = colSums(Biostrings::trinucleotideFrequency(hg38_seq))
  spectrumW = spectrum[names(spectrum) %in% SBS96_Subtypes]
  spectrumC = spectrum[!names(spectrum) %in% SBS96_Subtypes]
  names(spectrumC) = as.character(Biostrings::reverse(Biostrings::complement( Biostrings::DNAStringSet(names(spectrumC)) )))
  all(names(spectrumC) %in% SBS96_Subtypes)
  spectrum = matrix(rep(0,length(SBS32_Subtypes)),dimnames = list(SBS32_Subtypes,transcript))
  spectrum[names(spectrumW),] = spectrumW
  spectrum[names(spectrumC),] = spectrum[names(spectrumC),] + spectrumC
  spectrum = spectrum[SBS96_Subtypes,]
  return(spectrum)
}

#' get_SBS96_spectrum_driver
#'
#' @param driverlist a table containing a list of drivers alterations, with columns chr, pos, and alt
#'
#' @return A vector of size 96 containing the SBS spectrum of the transcript
#' @export
#'
#' @examples
#' TP53_LUAD.driver.spectrum = get_SBS96_driver_spectrum(TP53_drivers_intogen_LUAD)
#' print(TP53_LUAD.driver.spectrum)
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom Biostrings getSeq
#' @importFrom Biostrings trinucleotideFrequency
#' @importFrom Biostrings complement
#' @importFrom Biostrings reverse
#' @importFrom Biostrings DNAStringSet
get_SBS96_driver_spectrum = function( driverlist ){#gene="TP53",cohort="LUAD",observed=T){
  #@param gene the gene alias
  #@param cohort the cohort from intogen associated with the gene (examples: BRCA, PRAD)
  #@param observed boolean value specifying whether to restrict to observed mutations, or to use in silico predicted driver mutations from boostDM (default: TRUE)
  driver.preds = driverlist #intogen.boostDM.drivers %>% filter(gene==.env$gene,cohort%in%.env$cohort)
  #if(observed) driver.preds = driver.preds %>% filter(frequency>0)
  if(nrow(driver.preds)==0) stop("Gene-cohort pair not found in intogen data")

  SBS32_Subtypes = c("ACA", "ACC", "ACG", "ACT",
                     "CCA", "CCC", "CCG", "CCT",
                     "GCA", "GCC", "GCG", "GCT",
                     "TCA", "TCC", "TCG", "TCT",
                     "ATA", "ATC", "ATG", "ATT",
                     "CTA", "CTC", "CTG", "CTT",
                     "GTA", "GTC", "GTG", "GTT",
                     "TTA", "TTC", "TTG", "TTT")
  SBS96_Subtypes = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T","C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T","G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T","T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",
                     "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T","C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T","G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T","T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",
                     "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T","C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T","G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T","T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                     "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T","C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T","G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T","T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
                     "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T","C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T","G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T","T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
                     "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T","C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T","G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T","T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
  # retrieve transcript sequence + neighboring nucleotides from DB
  hg38_seq <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                                 names=driver.preds$chr ,
                                 start=driver.preds$pos-1,end=driver.preds$pos+1)) #,
                                 #strand=rtracklayer::strand(hg38_tr))
  hg38_seq.alt = driver.preds$alt
  # complement mutations
  to_comp = which(!hg38_seq %in% SBS32_Subtypes)
  hg38_seq[to_comp] = as.character(Biostrings::reverse(Biostrings::complement( Biostrings::DNAStringSet(hg38_seq[to_comp]) )))
  hg38_seq.alt[to_comp] = as.character(Biostrings::complement( Biostrings::DNAStringSet(hg38_seq.alt[to_comp]) ))
  # compute spectrum
  hg38_seq.mut = factor(paste0(substr(hg38_seq,1,1),"[",substr(hg38_seq,2,2),">",hg38_seq.alt,"]",substr(hg38_seq,3,3)),levels=SBS96_Subtypes)
  spectrum = table( hg38_seq.mut )

  return(spectrum)
}
