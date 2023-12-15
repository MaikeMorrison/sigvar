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
