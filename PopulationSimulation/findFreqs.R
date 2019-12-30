findFreqs <- function(pop){
  
  allelefreqs <- apply(pop,
                       MARGIN = 1, #rows are loci
                       FUN = mean) #amount of 1s/total # loci
  
  names(allelefreqs) <- paste(rep("locus", length(allelefreqs)), 1:length(allelefreqs), sep = ".")
  return(allelefreqs)
}

