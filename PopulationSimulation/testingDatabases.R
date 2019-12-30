library(dplyr)
library(DBI)
source("findcRMP.R")
source("popSimulation.R")
source("recombineFunction.R")
source("findFreqs.R")

conn <- dbConnect(RSQLite::SQLite(), "test.db")
dbListTables(conn)

  genome <- hypredGenome(num.chr = 1, len.chr = 4, num.snp.chr = 100)

population <- popSimulation(10,10, 10^-9, genome)
dbRemoveTable(conn, paste("test_population", 1,sep = "_"))
dbRemoveTable(conn, paste("allele_frequency", 1, sep = "_"))
dbWriteTable(conn, paste("test_population", 1,sep = "_"), as.data.frame(t(population[[1]]))) #replace the nummber with loop counter later 
dbWriteTable(conn, paste("allele_frequency", 1, sep = "_"), population[[2]]) #ditto ^
dbListTables(conn)


n <- sample(1:nrow(population[[1]]), 1) #random person from the population
#preallocate empty matrix for random individuals genome (each homolog is 1 row)
indiv <- matrix(nrow = 2, ncol = 100) # indiv <- matrix(nrow = 2, ncol = NUMBER OF SNPS)

#gets the homologs from the database
if(n %% 2 == 0){
#if the random number was even, then the 1st homolog is the previous one: previous one then RN  
  indiv[1, ] <- as.vector(t(dbGetQuery(conn, paste("SELECT V", n-1," FROM test_population_", 1, sep = ""))))
  indiv[2, ] <- as.vector(t(dbGetQuery(conn, paste("SELECT V", n," FROM test_population_", 1, sep = ""))))
} else {
#if RN was odd, RN is first homolog: RN then next one  
  indiv[1, ] <- as.vector(t(dbGetQuery(conn, paste("SELECT V", n," FROM test_population_", 1, sep = ""))))
  indiv[2, ] <- as.vector(t(dbGetQuery(conn, paste("SELECT V", n+1," FROM test_population_", 1, sep = ""))))
}

#gets last generation allele frequencies
popFreqs <- as.vector(t(dbGetQuery(conn, paste("SELECT Generation_", ncol(population[[2]])," FROM allele_frequency_", 1, sep = ""))))

findcRMP(popFreqs, indiv)
