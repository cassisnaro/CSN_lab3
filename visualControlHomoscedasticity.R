numberbins <- 7;

source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

#setEPS()#Set graphics as EPS

for(x in 1:nrow(source)){
  dataRead = read.table(source$file[x],header=FALSE)
  colnames(dataRead) = c("vertices","degree_2nd_moment","mean_length")
  language <- source$language[x]
  

  vars <- array(dim=c(numberbins))
  means <- array(dim=c(numberbins))
  tmpCut <- cut(dataRead$vertices, numberbins)
  dataRead$group <- as.numeric(tmpCut)
  mylevels <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", levels(tmpCut)) ),
        upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(tmpCut)) ))
  
  for (i in 1:numberbins){
    dataToProcess <- dataRead[which (dataRead$group==i),]
    means[i] <- mean(dataToProcess$vertices)
    vars[i] <- var(dataToProcess$degree_2nd_moment)
  }
  dataFrameCorrelation <- data.frame( Means = means, Vars = vars)
  cat(language,"\'s correlation: ",cor(dataFrameCorrelation, use ="complete.obs")["Means","Vars"],"");

  #postscript(paste('./homoscedasticity/',language,"_noModel",'.ps',sep = ""))
  plot(range(dataRead$vertices), range(0:max(vars, na.rm = TRUE)),type="n")
  for (i in 1:numberbins){
    segments(mylevels[i,1],vars[i],mylevels[i,2],vars[i])
    points(mean(mylevels[i,1],mylevels[i,2]),nrow(dataRead[which (dataRead$group==i),])/nrow(dataRead))
  }
  #dev.off()
  cat ("Press [enter] to continue")
  line <- readline()
  

#   if (x==10){
#     new_vars <- array(dim=c(numberbins))
#     for (i in 1:numberbins){
#       predicted_value <- (means[i]/2)^0.6736349;
#       new_vars[i] <- var(dataToProcess$degree_2nd_moment - predicted_value);
#     }
#     postscript(paste('./homoscedasticity/',language,"_model1",'.ps',sep = ""))
#     plot(range(dataRead$vertices), range(vars),type="n")
#     for (i in 1:numberbins){
#       dataToProcess <- dataRead[which (dataRead$group==i),]
#       segments(min(dataToProcess$vertices),new_vars[i],max(dataToProcess$vertices,new_vars[i]))
#     }  
#     dev.off()
#   }
}

