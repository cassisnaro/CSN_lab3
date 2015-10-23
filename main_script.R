tolerance <- 10E-6;
check_validity <- function(file) {
  dataRead = read.table(file,header=FALSE)
  numberRows <- dim(dataRead)[1]
  #Checking validity:
  validity <- TRUE;
  for( i in 1:numberRows){
    n = dataRead[i,1];
    meanKsquare = dataRead[i,2];
    meanD = dataRead[i,3];
    if(!((4-6/n)<=(meanKsquare + tolerance) && (meanKsquare-tolerance)<=(n-1) )){
      validity <- FALSE;
      break;
    }
    meanDLowerBound = n/(8*(n-1))*meanKsquare + 1/2;
    meanDUpperBound = n-1;
    if(!(meanDLowerBound<=(meanD + tolerance) && (meanD-tolerance)<=meanDUpperBound )){
      validity <- FALSE;
      break;
    }
  }
  if(validity){
    cat("The data in file ",file," is valid.\n");
  }else{
    cat("The data in file ",file," is unvalid.\n");
  }
}
write_summary <- function(language,file) {
  dataRead = read.table(file,header=FALSE);
  ns <- dataRead$V1;
  meanKsquares <- dataRead$V2;
  meanDs <- dataRead$V3;
  cat(language, dim(dataRead)[1], mean(ns), sd(ns), mean(meanKsquares), sd(meanKsquares),"\n");
}
preliminary_visualization <- function(language,file){
  languageData = read.table(file, header = FALSE);
  colnames(languageData) = c("vertices","degree_2nd_moment","mean_length")
  languageData = languageData[order(languageData$vertices), ]
  
  print(language)
  
  postscript(paste('./figures/',language,"_vertices","_meanLength",'.ps',sep = ""))
  plot(languageData$vertices, languageData$mean_length, xlab = "vertices", ylab = "mean dependency length", main = language)
  dev.off()
  
  postscript(paste('./figures/',language,"_logVertices","_logMeanLength",'.ps',sep = ""))
  plot(log(languageData$vertices), log(languageData$mean_length), xlab = "log(vertices)", ylab = "log(mean dependency length)", main = language)
  dev.off()
  
  mean_Language = aggregate(languageData, list(languageData$vertices), mean)
  postscript(paste('./figures/',language,"_meanVertices","_meanMeanLength",'.ps',sep = ""))
  plot(mean_Language$vertices, mean_Language$mean_length, xlab = "vertices", ylab = "mean mean dependency length", main = language)
  dev.off()
  
  postscript(paste('./figures/',language,"_logMeanVertices","_logMeanMeanLength",'.ps',sep = ""))
  plot(log(mean_Language$vertices), log(mean_Language$mean_length), xlab = "log(vertices)", ylab = "log(mean mean dependency length)", main = language)
  dev.off()
  
  postscript(paste('./figures/',language,"_logMeanVertices","_logMeanLength",'_plusEstimation','.ps',sep = ""))
  plot(log(mean_Language$vertices), log(mean_Language$mean_length), xlab = "log(vertices)", ylab = "log(mean mean dependency length)", main = language)
  lines(log(mean_Language$vertices),log(mean_Language$mean_length),col = "green")
  #lines(log(mean_Language$vertices),log((mean_Language$mean_vertices+1)/3),col = "red")
  
  dev.off()
}
source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)
for (x in 1:nrow(source)) {
  check_validity(source$file[x])
}
for (x in 1:nrow(source)) {
  write_summary(source$language[x], source$file[x])
}
for (x in 1:nrow(source)) {
  preliminary_visualization(source$language[x], source$file[x])
}