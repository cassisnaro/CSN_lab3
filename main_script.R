library(xtable)

#Initialization
source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)


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

#Check Validity
for (x in 1:nrow(source)) {
  check_validity(source$file[x])
}

#Summary
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


calcS <- function(m) {
  sqrt(deviance(m)/df.residual(m))
}

#For each language
#LanguageData <- Input data to fit.
ensemble_fitting <- function(languageData, language) {
  print(language)
  #AIC <-AIC(m)
  #s <- sqrt(deviance(m)/df.residual(m))
  #params <- coef(m)
  nModels <- 3
  
  AICv <- rep(-1,nModels)
  sv <- rep(-1,nModels)
  
  ##NonLinear Modelling
  #Model 1
  i <- 1
  b_initial_model_1 = log(mean(languageData$degree_2nd_moment)) / (log(mean(languageData$vertices))-log(2))
  cat("Initial values of the model 1: ", b_initial_model_1, "\n")  
  nonlinear_model_1 = nls(degree_2nd_moment~(vertices/2)^b,data=languageData, start = list( b = b_initial_model_1), trace = TRUE)
  
  AICv[i] <- AIC(nonlinear_model_1)
  sv[i] <- calcS(nonlinear_model_1)
  
  #Model 2
  i <- 2
  linear_model = lm(log(degree_2nd_moment)~log(vertices), languageData)
  a_initial_model_2 = coef(linear_model)[1]
  b_initial_model_2 = coef(linear_model)[2]
  cat("Initial values of the model 2: ", a_initial_model_2, " - ", b_initial_model_2, "\n")  
  nonlinear_model_2 = nls(degree_2nd_moment~a*vertices^b,data=languageData, start = list(a = a_initial_model_2, b = b_initial_model_2), trace = TRUE)
  
  AICv[i] <- AIC(nonlinear_model_2)
  sv[i] <- calcS(nonlinear_model_2)
  
  #Model 3
  i <- 3
  linear_model = lm(log(degree_2nd_moment)~vertices, languageData)
  a_initial_model_3 = exp(coef(linear_model)[1])
  c_initial_model_3 = coef(linear_model)[2]
  cat("Initial values of the model 3: ", a_initial_model_3, " - ", c_initial_model_3, "\n")  
  nonlinear_model_3 = nls(degree_2nd_moment~a*exp(vertices*c),data=languageData, start = list(a = a_initial_model_3, c = c_initial_model_3), trace = TRUE)
  
  AICv[i] <- AIC(nonlinear_model_3)
  sv[i] <- calcS(nonlinear_model_3)
  
  ##Plots
  #Model 1 plot
  postscript(paste('./figures/',language,"_model1",'.ps',sep = ""))
  plot(languageData$vertices, languageData$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment", main = language)
  b_model_1 = coef(nonlinear_model_1)["b"];
  lines(languageData$vertices, (languageData$vertices/2)^b_model_1, col="green")
  dev.off()
  
  #Model 2 plot
  postscript(paste('./figures/',language,"_model2",'.ps',sep = ""))
  plot(languageData$vertices, languageData$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment", main = language)
  a_model_2 = coef(nonlinear_model_2)["a"];
  b_model_2 = coef(nonlinear_model_2)["b"];
  lines(languageData$vertices, a_model_2*languageData$vertices^b_model_2, col="green")
  dev.off()
  
  #Model 3 plot
  postscript(paste('./figures/',language,"_model3",'.ps',sep = ""))
  plot(languageData$vertices, languageData$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment", main = language)
  a_model_3 = coef(nonlinear_model_3)["a"];
  c_model_3 = coef(nonlinear_model_3)["c"];
  lines(languageData$vertices, a_model_3*exp(languageData$vertices*c_model_3), col="green")
  lines(languageData$vertices, a_initial_model_3*exp(languageData$vertices*c_initial_model_3), col="red")
  dev.off()

  
  deltaAICv <- (AICv - AICv[which.min(AICv)])
  
  paramList <- list(
                    c(b_model_1),
                    c(a_model_2, b_model_2),
                    c(a_model_3, c_model_3)
                    )
  
  list(sv, AICv, deltaAICv, paramList)
}


#For each different language, read and apply
langData <- lapply(source$file, read_single_file)
modelResult <- mapply(function(x,y) ensemble_fitting(as.data.frame(x), y), langData, source$language)



#Create the 3 Table2 and Table 3
table2s <- data.frame()
table2AIC <- data.frame()
table2DeltaAIC <- data.frame()
table3 <- data.frame()
for (i in 1:length(source$language)) { #For each language
  #Warning, Gypsy programming ahead when creating the entries for each data.frame
  #x[[1]] <- s
  #x[[2]] <- AIC
  #x[[3]] <- delta AIC
  x <- modelResult[,i]
  table2s <- rbind(table2s, data.frame(Mod1=x[[1]][1], Mod2=x[[1]][2], Mod3=x[[1]][3]))
  table2AIC <- rbind(table2AIC, data.frame(Mod1=x[[2]][1], Mod2=x[[2]][2], Mod3=x[[2]][3]))
  table2DeltaAIC <- rbind(table2DeltaAIC, data.frame(Mod1=x[[3]][1], Mod2=x[[3]][2], Mod3=x[[3]][3]))
  
  #Parameter table
  #Model
  #Parameter (nested in model)
  par <- x[[4]]
  table3 <- rbind(table3, data.frame(
    b1= par[[1]][1], #Model 1
    a2=par[[2]][1], b2=par[[2]][2], #Model 2
    a3=par[[3]][1], c3=par[[3]][2] #Model 3
    ))
}

#change row names by language
rownames(table2s) <- source$language
rownames(table2AIC) <- source$language
rownames(table2DeltaAIC) <- source$language
rownames(table3) <- source$language




#Check homocedasticity!?

#source = read.table("list.txt", 
#                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
#                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
#)
#Generate models
#for (x in 1:nrow(source)) {
#  m <- model_fitting(source$language[x],source$file[x])
  #Check homocedasticity!?
  
#  RSS <- deviance(m)
#   AIC <- AIC(m)
#}

