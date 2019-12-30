library(hypred)


#workhorse function - recombines gametes
recombine <- function(genomeref, gamete1, gamete2, mutationrate){
  recombinedgamete <- hypredRecombine(genomeref,
                                      genomeA = gamete1,
                                      genomeB = gamete2,
                                      mutate = TRUE,
                                      mutation.rate.snp = mutationrate,
                                      mutation.rate.qtl = mutationrate,
                                      block = FALSE)
  return(recombinedgamete)
}


simulatepop <- function(G, N, mutRate, genome.ref) {
  founders <- hypredFounder(genome.ref, 1)
  
  loci <- genome.ref@num.chr * genome.ref@num.snp.chr
  #2 because diploid
  random.mate.pop <- matrix(nrow = N * 2, ncol = loci)
  
  #starter generation - can only get 1 from one parent and 0 from other so will be 1/0 at all loci
  random.mate.pop[seq(1, N*2, by = 2), ] <- founders[1, ]
  random.mate.pop[seq(2, N*2, by = 2), ] <- founders[2, ]
  
  # 4 because each diploid cell creates 4 haploid daughter cells via meiosis
  random.mate.pop.temp <- matrix(nrow = N * 2, ncol = loci)
  
  #pdf(paste("Plots/improved_allele_frequcy_trial_", trialnum, ".pdf", sep = ""))
  plot(
    x = 1:genome.ref@num.snp.chr,
    y = apply(random.mate.pop, MAR = 2, FUN = mean),
    type = "l",
    col = rgb(0, 0, 0, alpha = 1.0),
    ylim = c(0, 1),
    xlab = "SNP Number",
    ylab = "Allele Frequency",
    main = paste("N = ", N, " and G = ", G, sep = ""))
  
  mtext(paste("(Mutation rate =", mutRate, 
              ", Cromosome length = ", genome.ref@len.chr, 
              " (in M))", sep = ""), 
        side = 3, 
        line = 0)
  
  abline(h = c(1, 0), lty = 3, col = "grey")
  
  for (generation in 1:G) {
    ## loop over generations
    
    for (indiv in 1:N) {
      ## loop over individuals
      ## in order to keep the population stable, each parent pair will have 2 offspring
      ## for this to work, each parent needs to produce 2 gametes
      
      ## gamete 1
      random.mate.pop.temp[2*indiv-1,] <-
        recombine(genome.ref,
                  random.mate.pop[2 * indiv - 1, ],
                  random.mate.pop[2 * indiv - 0,],
                  mutRate)
      ## gamete 2
      random.mate.pop.temp[(2*indiv-0),] <-
        recombine(genome.ref, 
                  random.mate.pop[2*indiv-1, ], 
                  random.mate.pop[2*indiv-0,], 
                  mutRate)
      
    } ## end for N - individuals
    
    ## permutate - bumpin' uglies
    ## replace is false so no offspring can get the same gamete for both halfs of chromosome
    ## each offspring "chooses" their 2 parent gametes - consistent with W-F models
    random.mate.pop <- random.mate.pop.temp[sample(1 : (N * 2), replace = FALSE), ]
    View(random.mate.pop)
    # at this point, we have a new population
    lines(
      1:genome.ref@num.snp.chr,
      y = apply(random.mate.pop, MAR = 2, FUN = mean),
      col = rgb(0, 0, 0, alpha = 0.1)
    )
  } ## end for G
  #print(random.mate.pop)
  return(random.mate.pop)
  #dev.off() 
  #random.mate.pop
}


genome <- hypredGenome(num.chr = 1, 
                       len.chr = 10,
                       num.snp.chr = 100)

#current.population <- matrix(nrow = N*2, ncol = genome@num.chr * genome@num.snp.chr)
#current.population <- 
print(simulatepop(100, 10000, 10^-9, genome))

