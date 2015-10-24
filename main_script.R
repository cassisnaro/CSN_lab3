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
#Read one language file
read_single_file <- function(file) {
  languageData = read.table(file, header = FALSE);
  colnames(languageData) = c("vertices","degree_2nd_moment","mean_length")
  languageData[order(languageData$vertices), ]
}


#For each model
#modelFormula <- Formula of the non-linear model we want to use for fitting
#linearFormula <- Formula of linear model we will use to estimate initial parameters
#LanguageData <- Input data to fit.
model_fitting <- function(modelFormula, linearFormula, languageData){
  #a_initial <- 1
  #b_initial <- 1
  
  #Calculate using the a linear model
  linear_model = lm(linearFormula, languageData)
  a_initial = exp(coef(linear_model)[1])
  b_initial = coef(linear_model)[2]

  cat("Initial values of the model: ", a_initial, " - ", b_initial, "\n")  
  nls(modelFormula,data=languageData, cpstart = list(a = a_initial, b = b_initial), trace = TRUE)
}

calcS <- function(m) {
  sqrt(deviance(m)/df.residual(m))
}

#For each language
#LanguageData <- Input data to fit.
ensemble_fitting <- function(languageData) {
  modelNames <- c("2","2")
  nModels <- length(modelNames)
  #AIC <-AIC(m)
  #s <- sqrt(deviance(m)/df.residual(m))
  #params <- coef(m)
  
  AICv <- rep(0,nModels)
  sv <- rep(0,nModels)
 
  #model2Copy
  i <- 1
  nlmF <- degree_2nd_moment ~ a*vertices^b
  lmF <- log(degree_2nd_moment)~log(vertices)
  m1 <- model_fitting(nlmF,lmF,languageData)
  AICv[i] <- AIC(m1)
  sv[i] <- calcS(m1)
  
  #model2
  i <- 2
  nlmF <- degree_2nd_moment ~ a*vertices^b
  lmF <- log(degree_2nd_moment)~log(vertices)
  m2 <- model_fitting(nlmF,lmF,languageData)
  AICv[i] <- AIC(m2)
  sv[i] <- calcS(m2)
  
  deltaAICv <- (AICv - AICv[which.min(AICv)])
  
  list(sv, AICv, deltaAICv)
}


language <- source$language[1]
file <- source$file[1]


#For each different language, read and apply
langData <- lapply(source$file, read_single_file)
modelResult<- lapply(langData, ensemble_fitting)

#Create the 3 Table2
table2s <- data.frame()
table2AIC <- data.frame()
table2DeltaAIC <- data.frame()
for (i in 1:length(modelResult)) {
  #Warning, Gypsy programming ahead when creating the entries for each data.frame
  #x[[1]] <- s
  #x[[2]] <- AIC
  #x[[3]] <- delta AIC
  x <- modelResult[[i]]
  table2s <- rbind(table2s, data.frame(Mod1=x[[1]][1], Mod2=x[[1]][2]))
  table2AIC <- rbind(table2AIC, data.frame(Mod1=x[[2]][1], Mod2=x[[2]][2]))
  table2DeltaAIC <- rbind(table2DeltaAIC, data.frame(Mod1=x[[3]][1], Mod2=x[[3]][2]))

}


sapply(langData,View)

#Check homocedasticity!?

RSS <- deviance(m)
AIC <- AIC(m)
