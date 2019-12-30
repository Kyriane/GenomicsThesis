library(hypred)

recombine <- function(genomeref, gamete1, gamete2, mutationrate){
  recombinedgamete <- hypredRecombine(genomeref,
                                      genomeA = gamete1,
                                      genomeB = gamete2,
                                      mutate = TRUE,
                                      mutation.rate.snp = mutationrate,
                                      mutation.rate.qtl = mutationrate,
                                      block = FALSE)
  return(as.vector(recombinedgamete))
}
