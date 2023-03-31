#' A function to train isolation forest anomaly detector for QC data
#'
#' @param df_input comma-separated (.csv), metric file. It should contain a "peptide" column and the metrics columns. It should also include "Annotations" for each run.
#' @param predictions isolation forest prediction results.
#' @export
#' @import caret ggplot2 MASS dplyr
#' @import h2o
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- MSstatsQC::DataProcess(S9Site54)
#' pred<-MSstatsQC.UML.trainR(df_input = sampleData)
#' MSstatsQC.UML.testR(df_input = sampleData, predictions = pred)

MSstatsQC.UML.testR<-function(df_input, predictions, threshold=NULL){
  
  #Removing cols with non numeric vals
  pca_input  <- select_if(df_input, is.numeric)
  #Removing cols with constant variance
  pca_input <- pca_input[ , which(apply(pca_input, 2, var) != 0)]
  
  # Getting Principal Component
  pca <- prcomp(pca_input, scale. = TRUE , center = TRUE)
  pcainfo <-  summary(pca)
  
  x <- as.vector(pca[["x"]][,1])
  y <- as.vector(pca[["x"]][,2])
  z <- as.vector(pca[["x"]][,3])
  anom <- as.vector(predictions$anomaly_score)
  
  pca_df <-  as.data.frame(cbind(x,y,z,anom))
  
  # Ploting 3d Plot
  pca_df <- pca_df %>% mutate(anom = if_else(anom > threshold, 'Anomaly', 'Normal') )
  
  fig <- plot_ly(pca_df, x = ~x, y = ~y, z = ~z, color = ~anom, colors = c('#BF382A', '#0C4B8E'))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = paste('PC1 ', pcainfo$importance[2]*100,"%")),
                                     yaxis = list(title = paste('PC2 ', pcainfo$importance[5]*100,"%")),
                                     zaxis = list(title = paste('PC3 ', pcainfo$importance[8]*100,"%"))
  ))
  
  df_input <-  data$df 
  df_input$Predicted_Label <- pca_df$anom
  
  fig
}
