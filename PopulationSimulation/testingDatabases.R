library(dplyr)
library(DBI)
library(RSQLite)
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
dbWriteTable(conn, paste("test_population", 1,sep = "_"), (population[[1]])) #replace the nummber with loop counter later 
dbWriteTable(conn, paste("allele_frequency", 1, sep = "_"), population[[2]]) #ditto ^
dbListTables(conn)


n <- sample(1:10, 1) #random person from the population
## Replace 10 with N when making loop -> 10 is population size

#gets the homologs from the database
if(n %% 2 == 0){
  #if the random number was even, then the 1st homolog is the previous one: previous one then RN
  indiv <- dbGetQuery(conn, paste("SELECT V", n-1,", V", n, " FROM test_population_", 1, sep = ""))
} else {
  #if RN was odd, RN is first homolog: RN then next one
  indiv <- dbGetQuery(conn, paste("SELECT V", n,", V", n+1, " FROM test_population_", 1, sep = ""))
}


#gets last generation allele frequencies
popFreqs <- as.vector(t(dbGetQuery(conn, paste("SELECT Generation_", 10," FROM allele_frequency_", 1, sep = ""))))
# replace 10 with G when converting to loop
findcRMP(popFreqs, indiv)
