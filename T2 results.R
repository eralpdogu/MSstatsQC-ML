library(qcc)

g1<-ggplot(SimData, aes(AcquiredTime, RT)) + 
  geom_point()+
  geom_smooth(method="loess")+
  ylab("Retention time (sec)")+
  xlab("Run ID")+
  ylim(650, 800)+
  geom_vline(xintercept=25, color="red")+
  facet_wrap(~Precursor,scales = "free", ncol=4)+
  theme_light()

g2<-ggplot(SimData, aes(AcquiredTime, TotalArea, color=Precursor)) + 
    geom_line(show.legend = FALSE)+
    ylab("Total peak area")+
    xlab("Run ID")+
    #ylim(600, 1000)+
    geom_vline(xintercept=25, color="red")+
    theme_light()
  
g3<-ggplot(SimData, aes(AcquiredTime, FWHM, color=Precursor)) + 
    geom_line(show.legend = FALSE)+
    #ylim(600, 1000)+
  ylab("Full width at half maximum (FWHM)")+
  xlab("Run ID")+
  geom_vline(xintercept=25, color="red")+
  theme_light()
  
g4<-ggplot(SimData, aes(AcquiredTime, MassAccu, color=Precursor))+
  geom_line()+
  #ylim(600, 1000)+
  ylab("Mass Accuracy")+
  xlab("Run ID")+
  geom_vline(xintercept=25, color="red")+
  theme_light()+
  guides(color=guide_legend("Peptide"))

grid.arrange(g1,g2,g3,g4, ncol=2)

#T2 control chart
T2<-list()
for (i in 1:nlevels(Data.set$peptide)){
X<-scale(Data.set[Data.set$peptide==levels(Data.set$peptide)[i],3:6])
Ti<-mqcc(X[1:25,], type = "T2", newdata = X[1:50,], pred.limits = FALSE)
T2[[i]]<-Ti$newstats
for (j in 1:50) {if (T2[[i]][j]>Ti$limits[2]) {T2[[i]][j]<-T2[[i]][j]}}
}
T2<-data.frame(matrix(unlist(T2), nrow=50),stringsAsFactors=FALSE)
colnames(T2)<-levels(SimData$Precursor)
T2<-cbind(RUN=1:50, T2)

Results_melt <- melt(T2,id.vars ="RUN")
ggplot(Results_melt, aes(RUN, variable)) + 
  geom_tile(aes(fill = value), colour = "white") +
  labs(x ="Time", y=NULL)+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_gradient(low = "white", high = "blue",name = "T2 values")+
  ggtitle(label = NULL)+
  theme(legend.position="bottom", panel.background = element_blank(),
        plot.background = element_blank(), plot.margin = unit(c(0.1,0,0,0), "cm"),
        axis.ticks.length = unit(0, "pt"))

