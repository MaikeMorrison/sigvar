#' import_SigProfiler
#'
#' @param folder a folder containing SigProfilerExtractor results
#'
#' @return A list of tibbles containing the SigProfiler activity matrices
#' @export
#'
#' @examples
#' # Import dummy SigProfilerExtractor output results
#' SPfolder <- system.file("extdata", "SP", package = "sigvar")
#' Qlist <- import_SigProfiler(SPfolder)
#' print(Qlist)
#'
#' # Plot imported signatures
#' plot_dots(Qlist[[1]], normalized = FALSE)
#' @importFrom readr read_tsv
import_SigProfiler <- function(folder = ".") {
  # find activity files within the folder
  input_files <- list.files(folder, pattern = "Activities[_refit]*.txt", recursive = TRUE, full.names = TRUE)
  # get name (de novo or COSMIC)
  input_files.names <- vapply(
    list.files(folder,
      pattern = "Activities[_refit]*.txt",
      recursive = TRUE,
      full.names = FALSE
    ),
    function(x) {
      tmp <- strsplit(x, "/")[[1]]
      tmp <- tmp[length(tmp) - 2]
    },
    FUN.VALUE = character(1)
  )
  # if multiple solutions, keep only suggested solutions instead of all NMF solutions
  if (length(grep(input_files, pattern = "Suggested", value = FALSE)) > 0) {
    input_files.tokeep <- grep(input_files, pattern = "Suggested", value = FALSE)
  } else {
    input_files.tokeep <- seq_len(length(input_files))
  }
  # read files
  Qlist <- lapply(input_files[input_files.tokeep], read_tsv)
  # rename list entries
  names(Qlist) <- input_files.names[input_files.tokeep]
  return(Qlist)
}

#' get_SBS96_spectrum
#'
#' @param transcript the ensembl ID of the transcript
#' @param organism the name of the organism associated with the transcript (species)
#' @param ref_genome BSgenome reference genome
#'
#' @return A 96 x 1 matrix containing the SBS spectrum of the transcript
#' @export
#'
#' @examples
#' # Run on
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
#' library(ref_genome, character.only = TRUE)
#' spectrum <- get_SBS96_spectrum(transcript = "ENST00000269305.9", ref_genome=ref_genome)
#' print(spectrum)
#' @importFrom GenomicFeatures exons
#' @importFrom GenomeInfoDb seqlevelsStyle
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
#' @importFrom TxDb.Mmusculus.UCSC.mm10.knownGene TxDb.Mmusculus.UCSC.mm10.knownGene
get_SBS96_spectrum <- function(transcript = "ENST00000269305.9", organism = "Homo sapiens",ref_genome ) {
  if (sub("_| ", "", tolower(organism)) == "homosapiens") txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  if (sub("_| ", "", tolower(organism)) == "musmusculus") txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  if (!exists("txdb")) stop("Organism not found. Valid answers are Homo sapiens or Mus musculus.")
  SBS32_Subtypes <- c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT",
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT"
  )
  SBS96_Subtypes <- c(
    "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT",
    "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT",
    "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT",
    "ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT",
    "ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT",
    "ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT"
  )
  # retrieve transcript sequence + neighboring nucleotides from DB
  ref_tr <- GenomicFeatures::exons(txdb, filter = list(tx_name = transcript))

  # Get the reference genome
  tryCatch(
    error = function(cnd) {
      stop("Please provide the name of a valid BSgenome object.", call. = FALSE)
    },
    {
      ref_genome <- BSgenome::getBSgenome(ref_genome)
    }
  )
  
  GenomeInfoDb::seqlevelsStyle(ref_tr) <- GenomeInfoDb::seqlevelsStyle(ref_genome)
  ref_seq <- Biostrings::getSeq(ref_genome,
      names = rtracklayer::chrom(ref_tr),
      start = rtracklayer::start(ref_tr) - 1, end = rtracklayer::end(ref_tr) + 1,
      strand = rtracklayer::strand(ref_tr)
  )

  # compute spectrum
  spectrum <- colSums(Biostrings::trinucleotideFrequency(ref_seq))
  spectrumW <- spectrum[names(spectrum) %in% SBS96_Subtypes]
  spectrumC <- spectrum[!names(spectrum) %in% SBS96_Subtypes]
  names(spectrumC) <- as.character(Biostrings::reverse(Biostrings::complement(Biostrings::DNAStringSet(names(spectrumC)))))
  all(names(spectrumC) %in% SBS96_Subtypes)
  spectrum <- matrix(rep(0, length(SBS32_Subtypes)), dimnames = list(SBS32_Subtypes, transcript))
  spectrum[names(spectrumW), ] <- spectrumW
  spectrum[names(spectrumC), ] <- spectrum[names(spectrumC), ] + spectrumC
  spectrum <- spectrum[SBS96_Subtypes, ]
  return(spectrum)
}

#' get_SBS96_driver_spectrum
#'
#' @param driverlist a table containing a list of drivers alterations, with columns chr, pos, and alt
#' @param ref_genome BSgenome reference genome
#'
#' @return A a named integer table (length 96) containing the SBS spectrum of the transcript
#' @export
#'
#' @examples
#' data(TP53_drivers_intogen_LUAD, package = "sigvar")
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
#' library(ref_genome, character.only = TRUE)
#' TP53_LUAD.driver.spectrum <- get_SBS96_driver_spectrum(TP53_drivers_intogen_LUAD,ref_genome)
#' print(TP53_LUAD.driver.spectrum)
#' @importFrom Biostrings getSeq
#' @importFrom Biostrings trinucleotideFrequency
#' @importFrom Biostrings complement
#' @importFrom Biostrings reverse
#' @importFrom Biostrings DNAStringSet
get_SBS96_driver_spectrum <- function(driverlist, ref_genome ) {
  # @param driverlist a dataframe with the list of driver positions (columns chr and pos required)
  # @param ref_genome BSgenome reference genome corresponding to the driverlist positions
  driver.preds <- driverlist
  if (nrow(driver.preds) == 0) stop("Gene-cohort pair not found in intogen data")

  SBS32_Subtypes <- c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT",
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT"
  )
  SBS96_Subtypes <- c(
    "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",
    "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",
    "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
    "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
    "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
    "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
  )
  
  # Get the reference genome
  tryCatch(
    error = function(cnd) {
      stop("Please provide the name of a valid BSgenome object.", call. = FALSE)
    },
    {
      ref_genome <- BSgenome::getBSgenome(ref_genome)
    }
  )
  
  # retrieve transcript sequence + neighboring nucleotides from DB
  ref_seq <- as.character(Biostrings::getSeq(ref_genome,
      names = driver.preds$chr,
      start = driver.preds$pos - 1, end = driver.preds$pos + 1
    ))
  ref_seq.alt <- driver.preds$alt
  # complement mutations
  to_comp <- which(!ref_seq %in% SBS32_Subtypes)
  ref_seq[to_comp] <- as.character(Biostrings::reverse(Biostrings::complement(Biostrings::DNAStringSet(ref_seq[to_comp]))))
  ref_seq.alt[to_comp] <- as.character(Biostrings::complement(Biostrings::DNAStringSet(ref_seq.alt[to_comp])))
  # compute spectrum
  ref_seq.mut <- factor(paste0(substr(ref_seq, 1, 1), "[", substr(ref_seq, 2, 2), ">", ref_seq.alt, "]", substr(ref_seq, 3, 3)), levels = SBS96_Subtypes)
  spectrum <- table(ref_seq.mut)

  return(spectrum)
}
