findcRMP <- function(freqs, person){
  
  ## expects random person to be a locix2 data frame (2 columns, loci rows)
  #   with allele expression encoded as 1 and no expression as 0
  ## expects freqs to be a vector of allele (1) frequency in the population
  
  randomPerson <- person+1 #in order to not be deviding by 0s
  rmp <- vector(length = nrow(randomPerson))
  ## alleles ~ Binomial at each loci
  ## defining allele expression (ie hetero or dominant) as 'sucess'
  
  for(z in 1:nrow(randomPerson)){
    if(randomPerson[z, 1] / randomPerson[z, 2] == 1){ ## both copies have same allele -
      # means they're homo 1/1 or 2/2 = 1
      if(randomPerson[z, 1] == 2){ #checks if homo dom (i.e. 2/2)
        rmp[z] <- dbinom(2, 2, prob = freqs[z]) #homo
      }else{ #if not homo dom, has to be homo rec b/c both copies have same allele
        rmp[z] <- dbinom(0, 2, prob = freqs[z])
      }
    }else{
      rmp[z] <- dbinom(1, 2, prob = freqs[z]) #hetero
    }
  } #end for
  
  return(prod(rmp))
}