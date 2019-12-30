library(ggplot2)
genome <- hypredGenome(num.chr = 1,
                       len.chr = 10,
                       num.snp.chr = 50)

n <- 100
g <- 100

plot(x = 1:g, y = rep(0,g), 
     ylim = c(0,1), 
     col = rgb(0, 0, 0, alpha = 0),
     xlab = "Generation",
     ylab = "Proportion of Fixed (or lost) Loci",
     main = paste("N:", n, ", Loci loss over time", sep = ""))

for(i in 1:10){
  population <- popSimulation(G = g, N = n, mutRate = 2.5*10^-9, genome.ref = genome) #expects list
  gone <- apply(as.data.frame(population[[2]]), MAR = 2, FUN = function(x) length(which(x == 1 | x == 0))/(nrow(population[[2]])))
  lines(x = 1:length(gone), y = gone, col = rgb(0, 0, 0, alpha = 0.5))

}



