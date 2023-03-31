#' A function to train isolation forest anomaly detector for QC data
#'
#' @param df_input comma-separated (.csv), metric file. It should contain a "peptide" column and the metrics columns. It should also include "Annotations" for each run.
#' @param cf isolation parameter. Default assumes 10% of isolation.
#' @export
#' @import solitude ggplot2 dplyr
#' @import 
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- MSstatsQC::DataProcess(S9Site54)
#' head(sampleData)
#' MSstatsQC.UML.trainR(df_input = sampleData)

MSstatsQC.UML.trainR <- function(df_input, cf=0.10){
  
  df_input <-  df_input %>% select_if(~ !any(is.na(.)))
  #Initialising IF model
  model <- solitude::isolationForest$new(
    num_trees = 100,
    sample_size = base::round(nrow(df_input)*0.10 + 2),
    #replace = T,
    seed = 12345
  )
  
  # fitting on input 
  input<-model$fit(df_input)
  
  # Extracting Predictons 
  predictions <-  model$predict(df_input)
  conf <-  1-cf
  
  #A scertaining a threshold via boostrap
  n = length(predictions$anomaly_score)/2
  B = 10000
  boot_result = rep(NA, B)
  for (i in 1:B) {
    boot.sample = sample(n, replace = TRUE)
    boot_result[i] = quantile(predictions$anomaly_score[boot.sample], cf)
  }
  thresh<-median(boot_result)
  print(thresh)
  
  return(predictions)
}
