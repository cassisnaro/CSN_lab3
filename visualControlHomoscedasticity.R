numberbins <- 5;

source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

setEPS()#Set graphics as EPS

for(x in 1:nrow(source)){
  dataRead = read.table(source$file[x],header=FALSE)
  colnames(dataRead) = c("vertices","degree_2nd_moment","mean_length")
  language <- source$language[x]
  

  vars <- array(dim=c(numberbins))
  means <- array(dim=c(numberbins))
  dataRead$group <- as.numeric(cut(dataRead$vertices, numberbins))
  for (i in 1:numberbins){
    dataToProcess <- dataRead[which (dataRead$group==i),]
    means[i] <- mean(dataToProcess$vertices)
    vars[i] <- var(dataToProcess$degree_2nd_moment)
  }
  postscript(paste('./homoscedasticity/',language,"_noModel",'.ps',sep = ""))
  plot(range(dataRead$vertices), range(vars),type="n")
  for (i in 1:numberbins){
    dataToProcess <- dataRead[which (dataRead$group==i),]
    segments(min(dataToProcess$vertices),vars[i],max(dataToProcess$vertices,vars[i]))
  }  
  dev.off()
  

  if (x==10){
    new_vars <- array(dim=c(numberbins))
    for (i in 1:numberbins){
      predicted_value <- (means[i]/2)^0.6736349;
      new_vars[i] <- var(dataToProcess$degree_2nd_moment - predicted_value);
    }
    postscript(paste('./homoscedasticity/',language,"_model1",'.ps',sep = ""))
    plot(range(dataRead$vertices), range(vars),type="n")
    for (i in 1:numberbins){
      dataToProcess <- dataRead[which (dataRead$group==i),]
      segments(min(dataToProcess$vertices),new_vars[i],max(dataToProcess$vertices,new_vars[i]))
    }  
    dev.off()
  }
}

