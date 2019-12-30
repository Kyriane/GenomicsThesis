findFreqs <- function(pop){
  
  allelefreqs <- apply(pop,
                       MARGIN = 2, #columns are loci
                       FUN = mean) #amount of 1s/total # loci
  
  names(allelefreqs) <- paste(rep("locus", length(allelefreqs)), 1:length(allelefreqs), sep = ".")
  return(allelefreqs)
}
