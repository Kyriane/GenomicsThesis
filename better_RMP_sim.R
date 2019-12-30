library(hypred)

set.seed(33)
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

popSimulation <- function(G, N, mutRate, genome.ref){
  
  founders <- hypredFounder(genome.ref, 1)
  
  loci <- genome.ref@num.chr * genome.ref@num.snp.chr
  #2 because diploid
  random.mate.pop <- matrix(nrow = N*2, ncol = loci)
  
  #starter generation - can only get 1 from one parent and 0 from other so will be 1/0 at all loci
  random.mate.pop[seq(1, N*2, by = 2), ] <- founders[1, ]
  random.mate.pop[seq(2, N*2, by = 2), ] <- founders[2, ]
  
  random.mate.pop.temp <- matrix(nrow = N*2, ncol = loci)
  
  #pdf(paste("Plots/improved_allele_frequcy_trial_", trialnum, ".pdf", sep = ""))
  
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
    
    # at this point, we have a new population
    
  } ## end for G
  
  return(random.mate.pop)
  
}

findFreqs <- function(pop){
  
  freqs <- t(data.frame(apply(pop, #t changes it so that the colums are the loci
                              MARGIN = 2, 
                              FUN = function(x) length(which(x == 1)))/nrow(pop), 
                        apply(pop, 
                              MARGIN = 2, 
                              FUN = function(x) length(which(x == 0)))/nrow(pop)))
  rownames(freqs) <- c("SNP (1) freq", "Reference (0) freq")
  colnames(freqs) <- paste(rep("locus", ncol(freqs)), 1:ncol(freqs), sep = " ")
  return(freqs)
}   

findcRMP <- function(freqs, randomPerson){
  rmp <- vector(length = ncol(randomPerson))
  # alleles ~ Binomial
  for(z in 1:ncol(randomPerson)){
    if(randomPerson[1, z] / randomPerson[2, z] == 1){ # both copies have ame allele -
      # means they're homo 1/1 or 2/2 = 1
      if(randomPerson[1, z] == 2){ #checks if homo dom (i.e. 2/2)
        rmp[z] <- dbinom(2, 2, prob = freq.table[1,z]) #homo
        # 
      }else{ #if not homo dom, has to be homo rec b/c both copies have same allele
        rmp[z] <- dbinom(0, 2, prob = freq.table[1, z])
      }
    }else{
      rmp[z] <- dbinom(1, 2, prob = freq.table[1,z])
    }
  } # end for
  return(prod(rmp))
}


genome <- hypredGenome(num.chr = 1, 
                       len.chr = 10,
                       num.snp.chr = 100)

n <- 100
g <- 100
mutationRate <- 10^-9





snps <- seq(10, 500, by = 10)
cRMP_diff_snp_sizes <- snps

set.seed(33)
for(size in snps){
  genome <- hypredGenome(num.chr = 1, 
                         len.chr = 10,
                         num.snp.chr = size)
  
  simulatedPop <- popSimulation(g, n, mutationRate, genome)
  freq.table <- findFreqs(simulatedPop)
  
  individual <- sample(1:n, 1)
  # we take the 2 rows of that individual and +1 to every element
  # so we don't end up dividing by 0 when figuring out zygosity
  personGenome <- simulatedPop[c(individual*2-1, individual*2), ] + 1 
  cRMP_diff_snp_sizes[size/10] <- findcRMP(freq.table, personGenome)
}

plot(x = snps, y = -log10(cRMP_diff_snp_sizes), 
     xlab = "# of SNPs", type = "l",
     ylab = "-log(cRMP)",
     main = paste("G:", g, ", N:", n, sep = ""))
mtext(paste("µ:", mutationRate,", chrom length:10M, set.seed(33)", sep = ""), 
      side = 3,
      line = 0)
abline(lm(-log10(cRMP_diff_snp_sizes) ~ snps), col = rgb(0.8,0,0, alpha = 0.5))
text(400, 75, 
     paste("r^2 = ", 
           round(summary(lm(-log10(cRMP_diff_snp_sizes) ~ snps))$adj.r.squared, 6), sep = ""))

## diff pop sizes

gen.def <- hypredGenome(num.chr = 1,
                        len.chr = 10,
                        num.snp.chr = 20)

n <- seq(50, 1000, by = 50)
cRMP_diff_pop_sizes <- n

set.seed(33)
for (size in n){
  simulatedPop <- popSimulation(g, size, mutationRate, gen.def)
  freq.table <- findFreqs(simulatedPop)
  
  individual <- sample(1:size, 1)
  # we take the 2 rows of that individual and +1 to every element
  # so we don't end up dividing by 0 when figuring out zygosity
  personGenome <- simulatedPop[c(individual*2-1, individual*2), ] + 1 
  cRMP_diff_pop_sizes[size/50] <- findcRMP(freq.table, personGenome)
}

plot(x = n, y = -log10(cRMP_diff_pop_sizes), 
     xlab = "N", type = "l",
     ylab = "-log(cRMP)",
     main = paste("G:", g, ", #SNPs: 20", sep = ""))
mtext(paste("µ:", mutationRate,", chrom length:10M, set.seed(33)", sep = ""), 
      side = 3,
      line = 0)





vec1 <- matrix(nrow = 10, ncol = 2)
vec1[,c(1,2)] <- 1:10
vec2 <- rep(c(1,2), 5)

test <- apply(vec1, MAR = 2, FUN = function(x) length(which(x == vec2))/length(vec1[,1]))




