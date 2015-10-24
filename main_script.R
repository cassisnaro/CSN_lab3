library(xtable)

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
  data.frame(language=language, N=dim(dataRead)[1], MeanN=mean(ns), St.DevN=sd(ns), MeanX=mean(meanKsquares), St.DevX=sd(meanKsquares));
}

preliminary_visualization <- function(language,file){
  languageData = read.table(file, header = FALSE);
  colnames(languageData) = c("vertices","degree_2nd_moment","mean_length")
  languageData = languageData[order(languageData$vertices), ]
  
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
  
  postscript(paste('./figures/',language,"_logVertices","_logMeanMeanLength",'.ps',sep = ""))
  plot(log(mean_Language$vertices), log(mean_Language$mean_length), xlab = "log(vertices)", ylab = "log(mean mean dependency length)", main = language)
  dev.off()
  
  postscript(paste('./figures/',language,"_logVertices","_logMeanMeanLength",'_plusEstimation','.ps',sep = ""))
  plot(log(languageData$vertices), log(languageData$mean_length), xlab = "log(vertices)", ylab = "log(mean mean dependency length)", main = language)
  lines(log(mean_Language$vertices),log(mean_Language$mean_length),col = "green")
  lines(log(mean_Language$vertices),log((mean_Language$vertices+1)/3),col = "red")
  dev.off()
  
  postscript(paste('./figures/',language,"_vertices","_degree_2nd_moment",'_plusEstimation','.ps',sep = ""))
  plot(languageData$vertices, languageData$degree_2n_moment, xlab = "vertices", ylab = "degree 2nd moment", main = language)
  lines(mean_Language$vertices,mean_Language$degree_2nd_moment,col = "green")
  lines(languageData$vertices,
        (1-1/languageData$vertices)*(5-6/languageData$vertices),col = "red")
  lines(languageData$vertices, 4-6/languageData$vertices,col = "blue")
  lines(languageData$vertices, languageData$vertices-1,col = "blue")
  dev.off()
}

source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

for (x in 1:nrow(source)) {
  check_validity(source$file[x])
}

l <- data.frame()
for (x in 1:nrow(source)) {
  l <- rbind(l,write_summary(source$language[x], source$file[x]))
}
xtable(l)


#Plot
for (x in 1:nrow(source)) {
  preliminary_visualization(source$language[x], source$file[x])
}





#Ensemble of models
model_fitting <- function(language,file){
  languageData = read.table(file, header = FALSE);
  colnames(languageData) = c("vertices","degree_2nd_moment","mean_length")
  languageData = languageData[order(languageData$vertices), ]

  #a_initial <- 1
  #b_initial <- 1
  


  b_initial_model_1 = log(mean(languageData$degree_2nd_moment)) / (log(mean(languageData$vertices))-log(2))
  cat("Initial values of the model: ", b_initial_model_1, "\n")  
  nonlinear_model_1 = nls(degree_2nd_moment~(vertices/2)^b,data=languageData, start = list( b = b_initial_model_1), trace = TRUE)
  plot(languageData$vertices, languageData$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment", main = language)
  b_model = coef(nonlinear_model_1)["b"];
  lines(languageData$vertices, (languageData$vertices/2)^b_model, col="green")
  
  linear_model = lm(log(degree_2nd_moment)~log(vertices), languageData)
  a_initial = exp(coef(linear_model)[1])
  b_initial = coef(linear_model)[2]
  cat("Initial values of the model: ", a_initial, " - ", b_initial, "\n")  
  nonlinear_model_2 = nls(degree_2nd_moment~a*vertices^b,data=languageData, start = list(a = a_initial, b = b_initial), trace = TRUE)
  plot(languageData$vertices, languageData$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment", main = language)
  a_model = coef(nonlinear_model_2)["a"];
  b_model = coef(nonlinear_model_2)["b"];
  lines(languageData$vertices, a_model*languageData$vertices^b_model, col="green")

  linear_model = lm(log(degree_2nd_moment)~vertices, languageData)
  a_initial = exp(coef(linear_model)[1])
  c_initial = coef(linear_model)[2]
  cat("Initial values of the model: ", a_initial, " - ", b_initial, "\n")  
  nonlinear_model_3 = nls(degree_2nd_moment~a*exp(vertices*c),data=languageData, start = list(a = a_initial, c = c_initial), trace = TRUE)
  plot(languageData$vertices, languageData$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment", main = language)
  a_model = coef(nonlinear_model_3)["a"];
  c_model = coef(nonlinear_model_3)["c"];
  lines(languageData$vertices, a_model*exp(languageData$vertices*c_model), col="green")
  lines(languageData$vertices, a_initial*exp(languageData$vertices*c_initial), col="red")
}

language <- source$language[1]
file <- source$file[1]
m <- model_fitting(language,file)
#Check homocedasticity!?

RSS <- deviance(m)
AIC <- AIC(m)
