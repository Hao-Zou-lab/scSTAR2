#' This code mainly implements the process of OPLS-DA (Partial Least Squares Discriminant Analysis), and performs model training, visualization and significance testing based on the specified parameters and data.
#'
#' @param data1  data1 is a set of spectra from a biological class.
#' @param data2  data2 is a set of spectra from a biological class.
#' @param variables "variables"is the variable indexs.
#'
#' @return model,Integrate all results into a structure model, including the model's Q2 value, P value, permutation test results and other information.
OPLSDA = function(data1, data2, variables){
  # run OPLS-DA and display results
  # Input:
  # "data1" and "data2" are groups of spectra from two biological classes
  # "variables" is the variable indexs

  # The parameters are set in the OPLSDA para.txt file
  # nrcv - number of folds cross-validation
  # nc - number of correlated variables in X
  # ncox - number of orthogonal variables in X
  # ncoy - number of orthogonal variables in Y
  # r2cutoff - r2 cutoff value to identify discriminatory variables
  # p-cutoff - p cutoff value to identify discriminatory variables
  # np - number of permutation tests
  extdir<-system.file("extdata",package="scSTAR2")
  dirOPLSpara<-paste(extdir,"/OPLSDA para.txt",sep="")

  fid = read.table(dirOPLSpara,sep = "\t",header = T,row.names = 1)
  nrcv =  fid[1,1]
  nc = fid[2,1]
  ncox = fid[3,1]
  ncoy = fid[4,1]
  pCutoff = fid[5,1]
  np = fid[6,1]
  pCutoff = pCutoff/length(variables)


  ###format data
  X = rbind(data1,data2)
  Y = matrix(0,nrow(X),2)
  Y[1:nrow(data1),1] = 1
  Y[(nrow(data1)+1):nrow(X),2] = 1

  ##
  model <- o2pls_m(X, Y, nc, ncox, ncoy, nrcv, pCutoff)# permutation test
  Q2 <- model$Q2Yhatcum
  model <- o2pls_m(X, Y, nc, ncox, ncoy, 0, pCutoff)# OPLSDA without CV to identify significant variables
  model$Q2Yhatcum <- Q2
  model_1 <- dispopls_m(model,X,Y,variables) #plot results

  ##permutation test

  PT_Q2 <- PT(model, X, Y, nc, ncox, ncoy, nrcv, np)$PTQ2
  model$P_Q2 <- 1 - pnorm(model$Q2Yhatcum, mean = mean(PT_Q2), sd = sd(PT_Q2))

  sign_buffer_1<- model_1$Sign_buffer
  buffer_1=model_1$buffer_1
  buffer_1=buffer_1*1.1
  if (length(which(sign_buffer_1 < 0)) > 0 && length(which(sign_buffer_1 > 0)) == 0 && !is.null(model$P_Q2)) {
   text_data <- data.frame(x = 1, y = buffer_1, label = paste("p value=", as.character(model$P_Q2)))
    #x11()
   p4 <- model_1$p4
  (p8 <-p4 +
     geom_text(data = text_data, aes(x = x, y = y, label = label),size = 4, hjust = 0, vjust = 0))
  print(p8)
  }else if (length(which(sign_buffer_1 > 0)) > 0 && length(which(sign_buffer_1 < 0)) == 0 && !is.null(model$P_Q2)){
    text_data1 <- data.frame(x = 1, y = buffer_1, label = paste("p value=", as.character(model$P_Q2)))
    #x11()
  p3 <- model_1$p3
  (p7<-p3 +
     geom_text(data = text_data1, aes(x = x, y = y, label = label),size = 4, hjust = 0, vjust = 0))
  print(p7)
  }else if(length(which(sign_buffer_1 < 0)) > 0 && length(which(sign_buffer_1 > 0))>0 && !is.null(model$P_Q2)){
    text_data2 <- data.frame(x = 1, y = buffer_1 , label = paste("p value=", as.character(model$P_Q2)))
    #x11()
  p2 <- model_1$p2
  (p6 <- p2 +
      geom_text(data = text_data2, aes(x = x, y = y, label = label),size = 4, hjust = 0, vjust = 0))
  print(p6)
  }else{
    text_data3 <- data.frame(x = 1, y = buffer_1, label = paste("p value=", as.character(model$P_Q2)))
    #x11()
  p5 <- model_1$p5
  (p9 <-p5 +
      geom_text(data = text_data3, aes(x = x, y = y, label = label),size = 4, hjust = 0, vjust = 0))
  print(p9)
  }

 return(model)
}


