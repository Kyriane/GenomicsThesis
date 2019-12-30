findcRMP <- function(freqs, person){
  # expects random person to be a 2x#loci matrix with allele expression
  #    encoded as 1 and no expression encoded as 0
  # expects freqs to be a vector of allele (1) freqency in the population
  randomPerson <- person+1 #in order to not be deviding by 0s
  rmp <- vector(length = ncol(randomPerson))
  # alleles ~ Binomial
  # defining allele expression (ie hetero or dominant) as 'sucess'
  for(z in 1:ncol(randomPerson)){
    if(randomPerson[1, z] / randomPerson[2, z] == 1){ # both copies have same allele -
      # means they're homo 1/1 or 2/2 = 1
      if(randomPerson[1, z] == 2){ #checks if homo dom (i.e. 2/2)
        rmp[z] <- dbinom(2, 2, prob = freqs[z]) #homo
        # 
      }else{ #if not homo dom, has to be homo rec b/c both copies have same allele
        rmp[z] <- dbinom(0, 2, prob = freqs[z])
      }
    }else{
      rmp[z] <- dbinom(1, 2, prob = freqs[z])
    }
  } # end for
  return(prod(rmp))
}