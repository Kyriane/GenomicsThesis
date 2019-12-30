library(dplyr)
source("findcRMP.R")
source("popSimulation.R")
source("parralellePopSim.R")
source("recombineFunction.R")
source("findFreqs.R")

genome2 <- hypredGenome(1,10,100)

#testing my functions as a whole
system.time({popSimulation(100, 1000, 10^-9, genome2)})

registerDoParallel(cores = 24)
system.time({popSimulationParallel(100, 1000, 10^-9, genome2)})



#need this data setup to test the loop
N <- 100
random.mate.pop <- matrix(nrow = 2*N, ncol = 20)
random.mate.pop[c(TRUE, FALSE), ] <- rep(1, 20)
random.mate.pop[c(FALSE, TRUE), ] <- rep(0, 20)

random.mate.pop.temp <- matrix(nrow = N*2, ncol = 20)
test <- matrix(nrow = N, ncol = 20)

#testing just the loop part of the function
system.time({
  test <-
    foreach(strand = seq(1, 2 * N, by = 2)) %dopar% recombine(genome2,
                                                              random.mate.pop[strand,],
                                                              random.mate.pop[strand +1,],
                                                              10^-9)
})


#testing non-parallel loop
system.time({
  for (indiv in 1:N) {
    random.mate.pop.temp[(2*indiv-0),] <-
      recombine(genome2, 
                random.mate.pop[2*indiv-1, ], 
                random.mate.pop[2*indiv-0,], 
                10^-9)
  } 
})


founders <- hypredFounder(genome2, 1)

loci <- genome2@num.chr * genome2@num.snp.chr
#2 because diploid
random.mate.pop <- matrix(nrow = N*2, ncol = loci)

#testing difference b/w 2 definition methods - spoiler alert: they're ~ the same
system.time({random.mate.pop[c(TRUE, FALSE), ] <- founders[1, ]})

system.time({random.mate.pop[seq(2, N*2, by = 2), ] <- founders[2, ]})




