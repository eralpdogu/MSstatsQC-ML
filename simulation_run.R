library(readxl)
library(h2o)
library(caret)
library(MASS)

factorial <- read_xlsx("Factorial Runs/Factorialcombinatins.xlsx",sheet = 1)

source("add_features.R")
source("ml_algo.R")


for (sim.size in c(10, 25, 50, 500)){
  num_peptides <- nlevels(guide.set.scale$peptide)
  for(i in 1:nrow(factorial)){
    sample_data<- list()
    
    ###### In control observation ~ 5* sim size  the of the actual 
    for(k in 1:nlevels(guide.set.scale$peptide)){ 
      sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k], sim.size*5)
      sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size*5),
                             sample_data_k,
                             RESPONSE= c("PASS"))
      sample_data[[k]] <- sample_data_k
    }
    data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
    data <- data.set[, !duplicated(colnames(data.set))]
    
    
    for(j in 2:6){
      #small mean shift up for all peptides all metrics
      if(factorial[i,j]== "-" & colnames(factorial[i,j])=="A"){  
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                                 sample_data_k[1]+0.5*IQR(sample_data_k[,1]), 
                                 sample_data_k[2]+0.5*IQR(sample_data_k[,2]), 
                                 sample_data_k[3]+0.5*IQR(sample_data_k[,3]), 
                                 sample_data_k[4]+0.5*IQR(sample_data_k[,4]),
                                 RESPONSE= c("FAIL"))
          sample_data[[k]] <- sample_data_k
          
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large mean shift up for all peptides all metrics
      else if(factorial[i,j]== "+" & colnames(factorial[i,j])=="A"){  
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                                 sample_data_k[1]+2*IQR(sample_data_k[,1]), 
                                 sample_data_k[2]+2*IQR(sample_data_k[,2]), 
                                 sample_data_k[3]+2*IQR(sample_data_k[,3]), 
                                 sample_data_k[4]+2*IQR(sample_data_k[,4]),
                                 RESPONSE= c("FAIL"))
          sample_data[[k]] <- sample_data_k
          
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large mean shift down for all peptides all metrics
      else if(factorial[i,j]== "+" & colnames(factorial[i,j])=="B"){  
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                                 sample_data_k[1]-2*IQR(sample_data_k[,1]), 
                                 sample_data_k[2]-2*IQR(sample_data_k[,2]), 
                                 sample_data_k[3]-2*IQR(sample_data_k[,3]), 
                                 sample_data_k[4]-2*IQR(sample_data_k[,4]),
                                 RESPONSE= c("FAIL"))
          sample_data[[k]] <- sample_data_k
          
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #small mean shift down for all peptides all metrics
      else if(factorial[i,j]== "-" & colnames(factorial[i,j])=="B"){  
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                                 sample_data_k[1]-0.5*IQR(sample_data_k[,1]), 
                                 sample_data_k[2]-0.5*IQR(sample_data_k[,2]), 
                                 sample_data_k[3]-0.5*IQR(sample_data_k[,3]), 
                                 sample_data_k[4]-0.5*IQR(sample_data_k[,4]),
                                 RESPONSE= c("FAIL"))
          sample_data[[k]] <- sample_data_k
          
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large mean drift up for all peptides all metrics
      else if(factorial[i,j]== "+" & colnames(factorial[i,j])=="C"){  
        sample <- list()
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[k],1),
                                                           sample_data[m,1]+log(m,base=10)*IQR(sample_data[,1]), 
                                                           sample_data[m,2]+log(m,base=10)*IQR(sample_data[,2]),
                                                           sample_data[m,3]+log(m,base=10)*IQR(sample_data[,3]),
                                                           sample_data[m,4]+log(m,base=10)*IQR(sample_data[,4]),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", "RT", "TotalArea", "MassAccu", "FWHM","RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #small mean drift up for all peptides all metrics
      else if(factorial[i,j]== "-" & colnames(factorial[i,j])=="C"){  
        sample <- list()
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[k],1),
                                                           sample_data[m,1]+log(m,base=2)*IQR(sample_data[,1]), 
                                                           sample_data[m,2]+log(m,base=2)*IQR(sample_data[,2]),
                                                           sample_data[m,3]+log(m,base=2)*IQR(sample_data[,3]),
                                                           sample_data[m,4]+log(m,base=2)*IQR(sample_data[,4]),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", "RT", "TotalArea", "MassAccu", "FWHM","RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large mean drift down for all peptides all metrics
      else if(factorial[i,j]== "+" & colnames(factorial[i,j])=="D"){  
        sample <- list()
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[k],1),
                                                           sample_data[m,1]-log(m,base=10)*IQR(sample_data[,1]), 
                                                           sample_data[m,2]-log(m,base=10)*IQR(sample_data[,2]),
                                                           sample_data[m,3]-log(m,base=10)*IQR(sample_data[,3]),
                                                           sample_data[m,4]-log(m,base=10)*IQR(sample_data[,4]),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", "RT", "TotalArea", "MassAccu", "FWHM","RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #small mean drift down for all peptides all metrics
      else if(factorial[i,j]== "-" & colnames(factorial[i,j])=="D"){  
        sample <- list()
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[k],1),
                                                           sample_data[m,1]-log(m,base=2)*IQR(sample_data[,1]), 
                                                           sample_data[m,2]-log(m,base=2)*IQR(sample_data[,2]),
                                                           sample_data[m,3]-log(m,base=2)*IQR(sample_data[,3]),
                                                           sample_data[m,4]-log(m,base=2)*IQR(sample_data[,4]),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", "RT", "TotalArea", "MassAccu", "FWHM","RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large Cyclic pattern
      else if(factorial[i,j]== "+" & colnames(factorial[i,j])=="E"){
        sample <- list()
        for(k in 1:nlevels(guide.set$peptide)){
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k], sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[j],1),
                                                           sample_data[m,1]+2*sin(2*pi*m/sim.size), 
                                                           sample_data[m,2]+2*sin(2*pi*m/sim.size), 
                                                           sample_data[m,3]+2*sin(2*pi*m/sim.size),
                                                           sample_data[m,4]+2*sin(2*pi*m/sim.size),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", "RT", "TotalArea", "MassAccu", "FWHM","RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      # small Cyclic pattern
      else if(factorial[i,j]== "-" & colnames(factorial[i,j])=="E"){
        sample <- list()
        for(k in 1:nlevels(guide.set$peptide)){
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k], sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[j],1),
                                                           sample_data[m,1]+sin(2*pi*m/sim.size), 
                                                           sample_data[m,2]+sin(2*pi*m/sim.size), 
                                                           sample_data[m,3]+sin(2*pi*m/sim.size),
                                                           sample_data[m,4]+sin(2*pi*m/sim.size),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", "RT", "TotalArea", "MassAccu", "FWHM","RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
    }# column ends 
    data <- data[sample(nrow(data), nrow(data)), ] # shuffle the data
    
    
    
    dl_model <- ml_algo(data,i)
    
    
    
    cf<- data.frame(h2o.confusionMatrix(dl_model,valid = T),stringsAsFactors = F)
    factorial[i,"Overall error rate"] <- (1-cf[3,3])*100
    factorial[i,"False positives"] <- cf[2,3]
    factorial[i,"False negatives"] <- cf[1,3]
    write.csv(factorial, file = paste("factorial_guideset_",nrow(guide.set.scale),"_",sim.size,".csv",sep = ""))
    
  }
  
  
  
  
  
}
