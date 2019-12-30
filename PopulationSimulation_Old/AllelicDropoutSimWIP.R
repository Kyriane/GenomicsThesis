library(hypred)
library(ggplot2)
## Test for Allelic Dropout/Fixation

genome.ref <- hypredGenome(num.chr = 1,
                           len.chr = 10,
                           num.snp.chr = 100)
G = 100
N = 50
mutRate <- 2.5*10^-9
founders <- hypredFounder(genome.ref, 1)

loci <- genome.ref@num.chr * genome.ref@num.snp.chr
#2 because diploid
random.mate.pop <- matrix(nrow = N*2, ncol = loci)

#starter generation - can only get 1 from one parent and 0 from other so will be 1/0 at all loci
random.mate.pop[seq(1, N*2, by = 2), ] <- founders[1, ]
random.mate.pop[seq(2, N*2, by = 2), ] <- founders[2, ]

random.mate.pop.temp <- matrix(nrow = N*2, ncol = loci)

#pdf(paste("Plots/improved_allele_frequcy_trial_", trialnum, ".pdf", sep = ""))
amountDropped <- data.frame(generation0 = rep(0.5, loci))

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

set.seed(1)
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
  amountDropped <- cbind(amountDropped, amount = apply(random.mate.pop, MAR = 2, FUN = mean))
  
  lines(
    1:genome.ref@num.snp.chr,
    y = apply(random.mate.pop, MAR = 2, FUN = mean),
    col = rgb(0, 0, 0, alpha = 0.1)
  )
  
} ## end for G
dropped <- as.data.frame(t(amountDropped))
colnames(dropped) <- paste(rep("SNP", ncol(dropped)), 1:ncol(dropped), sep = ".")
rownames(dropped) <- paste(rep("Generation", nrow(dropped)), 0:(nrow(dropped)-1), sep = ".")
colour = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] #makes vec w colours - greys

plot(1:nrow(dropped), dropped$SNP.1, type = "l", ylim = c(0,1), xlab = "Generation", ylab = "frequency")

for(i in 2:nrow(dropped)){
lines(x = 1:nrow(dropped), y = dropped[, i], col = sample(colour, 1))
}
for(i in 1:nrow(dropped)){
  print(length(which(dropped[i,] < 1 & dropped[i, ] > 0)))
}

gone <- apply(dropped, MAR = 1, FUN = function(x) length(which(x < 1 & x > 0))/nrow(dropped))
plot(1:length(gone), gone, type = "l", ylim = c(0,1))
