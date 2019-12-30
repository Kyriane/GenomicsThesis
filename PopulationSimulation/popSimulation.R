library(hypred)

popSimulation <- function(G, N, mutRate, genome.ref){
  
  founders <- hypredFounder(genome.ref, 1)
  
  
  loci <- genome.ref@num.chr * genome.ref@num.snp.chr
  #2 because diploid
  #every 2 columns represents a person -> 2 homologs
  random.mate.pop <- as.data.frame(matrix(nrow = loci, ncol = N*2))
  
  #starter generation - can only get 1 from one parent and 0 from other so will be 1/0 at all loci
  random.mate.pop[ , seq(1, N*2, by = 2)] <- founders[1, ]
  random.mate.pop[ , seq(2, N*2, by = 2)] <- founders[2, ]
  
  random.mate.pop.temp <- random.mate.pop
  
  propExpressedPerGen <- as.data.frame(matrix(nrow = loci, ncol = G)) #holds 1 frequencies for each gen
  
  colnames(propExpressedPerGen) <- paste(rep("Generation", ncol(propExpressedPerGen)), 
                                         1:(ncol(propExpressedPerGen)), 
                                         sep = "_") #names the columns

  #pdf(paste("Plots/improved_allele_frequcy_trial_", trialnum, ".pdf", sep = ""))
  
  for (generation in 1:G) {
    ## loop over generations
    
    for (indiv in 1:N) {
      ## loop over individuals
      ## in order to keep the population stable, each parent pair will have 2 offspring
      ## for this to work, each parent needs to produce 2 gametes
      
      ## gamete 1
      random.mate.pop.temp[ ,(2*indiv-1)] <-
        recombine(genome.ref,
                  random.mate.pop[ ,(2 * indiv - 1)],
                  random.mate.pop[ ,(2 * indiv - 0)],
                  mutRate)
      ## gamete 2
      random.mate.pop.temp[ ,(2*indiv-0)] <-
        recombine(genome.ref, 
                  random.mate.pop[ ,(2*indiv-1)] , 
                  random.mate.pop[ ,(2*indiv-0)], 
                  mutRate)
      
    } ## end for N - individuals
    
    ## permutate - bumpin' uglies
    ## replace is false so no offspring can get the same gamete for both halfs of chromosome
    ## each offspring "chooses" their 2 parent gametes - consistent with W-F models
    random.mate.pop <- random.mate.pop.temp[ ,sample(1 : (N * 2), replace = FALSE)]
    
    # at this point, we have a new population
    propExpressedPerGen[ , generation] <- findFreqs(random.mate.pop)
    
  } ## end for G
  
  propExpressedPerGen <- as.data.frame(Generation_0 = rep(0.5, loci), propExpressedPerGen) #adds gen0 (ie first heterozygotes)
  
  values <- list(random.mate.pop, propExpressedPerGen)
  return(values)
  
}