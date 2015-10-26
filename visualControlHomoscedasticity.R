numberbins <- 50;

source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

for(x in 1:nrow(source)){
  dataRead = read.table(source$file[x],header=FALSE)
  language <- source$language[x]
  

  vars <- array(dim=c(numberbins))
  means <- array(dim=c(numberbins))
  dataRead$group <- as.numeric(cut(dataRead$V1, numberbins))
  postscript(paste('./homoscedasticity/',language,"_noModel",'.ps',sep = ""))
  plot(range(dataRead$V1), range(0,5),type="n")
  for (i in 1:numberbins){
    dataToProcess <- dataRead[which (dataRead$group==i),]
    means[i] <- mean(dataToProcess$V1)
    vars[i] <- var(dataToProcess$V2)
    segments(min(dataToProcess$V1),vars[i],max(dataToProcess$V1,vars[i]))
  }
  dev.off()
}

