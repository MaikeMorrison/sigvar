import_SigProfiler <- function(folder="."){
  # find activity files within the folder
  input_files        = list.files(folder,pattern = "Activities[_refit]*.txt",recursive = T,full.names = T)
  # get name (de novo or COSMIC)
  input_files.names  = sapply( list.files(folder,pattern = "Activities[_refit]*.txt",recursive = T,full.names = F), function(x){strsplit(x ,"/")[[1]][3]} )
  # keep only suggested solutions instead of all NMF solutions
  input_files.tokeep = grep(input_files, pattern = "Suggested" ,value = F)
  # read files
  Qlist = lapply( input_files[input_files.tokeep], read.table,h=T)
  # rename list entries
  names(Qlist) = input_files.names[input_files.tokeep]
  return(Qlist)
}
