library(hypred)
library(parallel)

popSimulationParallel <- function(G, N, mutRate, genome.ref){
  
  founders <- hypredFounder(genome.ref, 1)
  
  loci <- genome.ref@num.chr * genome.ref@num.snp.chr
  #2 because diploid
  random.mate.pop <- matrix(nrow = N*2, ncol = loci)
  
  #starter generation - can only get 1 from one parent and 0 from other so will be 1/0 at all loci
  random.mate.pop[c(TRUE, FALSE), ] <- founders[1, ]
  random.mate.pop[c(FALSE, TRUE), ] <- founders[2, ]
  
#  random.mate.pop.temp <- matrix(nrow = N*2, ncol = loci)
  
  propExpressedPerGen <- as.data.frame(matrix(nrow = loci, ncol = G)) #hold 1 frequencies for each gen
  colnames(propExpressedPerGen) <- paste(rep("Generation", ncol(propExpressedPerGen)), 
                                         1:(ncol(propExpressedPerGen)), 
                                         sep = ".")  #names the columns
  
  #pdf(paste("Plots/improved_allele_frequcy_trial_", trialnum, ".pdf", sep = ""))
  
  for (generation in 1:G) {
    ## loop over generations
    
    ## loop over individuals
    ## in order to keep the population stable, each parent pair will have 2 offspring
    ## for this to work, each parent needs to produce 2 gametes
    
    # random.mate.pop.temp[seq(1, 2 * N, 2),] <-
    #   foreach(strand = seq(1, 2 * N, by = 2),
    #           .combine = rbind) %dopar% recombine(genome.ref,
    #                                               random.mate.pop[strand,],
    #                                               random.mate.pop[strand +1,],
    #                                               mutRate)
    # random.mate.pop.temp[seq(2, 2 * N, 2),] <-
    #   foreach(strand = seq(2, 2 * N, by = 2),
    #           .combine = rbind) %dopar% recombine(genome.ref,
    #                                               random.mate.pop[strand,],
    #                                               random.mate.pop[strand -1,],
    #                                               mutRate)
    # 
    
    # going to use mclapply() instead - loop over individuals, 1:N
    res <- mclapply(1:N, mc.cores = 32, FUN = function(y) {
      gamete1 <- recombine(genome.ref,
                           random.mate.pop[2 * y - 1, ],
                           random.mate.pop[2 * y - 0,],
                           mutRate)
      gamete2 <- recombine(genome.ref, 
                           random.mate.pop[2 * y - 1, ], 
                           random.mate.pop[2* y - 0,], 
                           mutRate)
      return(rbind(gamete1, gamete2))
    })
    #require(dplyr)
    #random.mate.pop.temp <- dplyr::bind_rows(res, .id = 'source')
    random.mate.pop.temp <- do.call("rbind", res)
    ## end for N - individuals
    
    ## permutate - bumpin' uglies
    ## replace is false so no offspring can get the same gamete for both halfs of chromosome
    ## each offspring "chooses" their 2 parent gametes - consistent with W-F models
    random.mate.pop <-
      random.mate.pop.temp[sample(1:(N * 2), replace = FALSE),]
    # at this point, we have a new population
    
    # grab the 1 allele expression for each loci at this generation
    propExpressedPerGen[ , generation] <- findFreqs(random.mate.pop) 
    
  } ## end for G
  
  propExpressedPerGen <- cbind(Generation.0 = rep(0.5, loci), propExpressedPerGen) #adds gen0 (ie first heterozygotes)
  
  return(list(random.mate.pop, propExpressedPerGen))
  
}
