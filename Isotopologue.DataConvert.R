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
#' head(sampleData)
#' # Find the name of the peptides
#' levels(sampleData$Precursor)
#' # Calculate change point statistics
#' MSstatsQC.UML.trainR(df_input = sampleData[1:20,])

Isotopologue.DataConvert<-function(Isotop.Data, transform=NULL){
  Isotop.Data.Final<-NULL
  for(i in 1:nlevels(Isotop.Data$peptide)) {
    Isotop.Data.Peptide<-Isotop.Data[Isotop.Data$peptide==levels(Isotop.Data$peptide)[i],]
    if (transform=="log") {
      model<-lm(log(Isotop.Data.Peptide$HeavyArea)~Isotop.Data.Peptide$Concentration)
      Isotop.Data.Peptide$HeavyArea<-scale(model$residuals)
    }
    else if (transform=="sqrt") {
      model<-lm(sqrt(Isotop.Data.Peptide$HeavyArea)~Isotop.Data.Peptide$Concentration)
      Isotop.Data.Peptide$HeavyArea<-scale(model$residuals)
    }
    else {
      model<-lm(Isotop.Data.Peptide$HeavyArea~Isotop.Data.Peptide$Concentration)
      Isotop.Data.Peptide$HeavyArea<-scale(model$residuals)
    }
    Isotop.Data.Peptide$RetentionTime<-rnorm(dim(Isotop.Data.Peptide)[1], 
                                             mean(Isotop.Data.Peptide$RetentionTime), 1)
    Isotop.Data.Final<-rbind(Isotop.Data.Final, Isotop.Data.Peptide)
  }
  Isotop.Data<-Isotop.Data.Final
  return(Isotop.Data)
}
