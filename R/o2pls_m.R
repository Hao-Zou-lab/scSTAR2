#'O2-PLS function used to build the model
#'For further predictions use the o2pls_pred function
#'
#' @param X matrix input.
#' @param Y matrix input.
#' @param nc number of correlated components, it can be determinated by PCA of Y'X.
#' @param ncox number of orthogonal components in X.
#' @param ncoy number of orthogonal component in Y.
#' @param nrcv number of fold in the cross-validation (full cross validation).
#' @param pCutoff p cutoff value to identify discriminatory variables.
#'
#' @return model,The return value contains various parameters and results of the fitted PLSR or O2PLS model, which can be used for subsequent analysis and evaluation of model performance.


o2pls_m =  function(X,Y,nc,ncox,ncoy,nrcv, pCutoff){
  #
  # O2-PLS function used to build the model
  # For further predictions use the o2pls_pred function
  #
  # Input:
  # X - X matrix input
  # Y - Y matrix input
  # prep - preprocessing available
  #           - 'no' : no preprocessing
  #           - 'mc' : meancentering
  #           - 'uv' : Univariance scaling
  #
  # nc - number of correlated components, it can be determinated by PCA of Y'X
  # ncox - number of orthogonal components in X
  # ncoy - number of orthogonal component in Y
  # nrcv - number of fold in the cross-validation (full cross validation)
  # v1.0, Xin Zou, Kent
  
  muClass1 = colMeans(X[which(Y[,1]==1),])
  muClass2 = colMeans(X[which(Y[,2]==1),])
  Signs = sign(muClass2-muClass1)

  nsx = nrow(X)
  nvx = ncol(X)
  nsy = nrow(Y)
  nvy = ncol(Y)

  model = list()
  model$ns = nsx

  model$nc=nc
  model$ncox=ncox
  model$ncoy=ncoy

  if(nsx!=nsy){
    print('Number of samples are different in X and Y')
    model = list()
  }

  model$MeanX = colMeans(X)
  model$MeanY = colMeans(Y)

  model$StdX=wrMisc::colSds(X)
  model$StdY=wrMisc::colSds(Y)

  model$SSX= sum(colSums(X^2))
  model$SSY=sum(colSums(Y^2))
  model$CSX=colSums(X^2)
  model$CSY=colSums(Y^2)

  if(nrcv == 0){
    #unite variance scaling
    X = X-t(matrix(model$MeanX, nrow = length(model$MeanX), ncol = nsx))
    Y = Y-t(matrix(model$MeanY, nrow = length(model$MeanY), ncol = nsy))
    X = X/t(matrix(model$StdX, nrow = length(model$StdX), ncol = nsx))
    Y = Y/t(matrix(model$StdY, nrow = length(model$StdY), ncol = nsy))

    M <- mjrO2pls(X,Y,nc,ncox,ncoy,'standard')
    M$Q2Yhatcum <- numeric(0)
    model <- M

    # pearson correlation between each variable X and PLS score model.T
    nsx <- nrow(X)
    nvx <- ncol(X)
    buffer <- matrix(rep(model$T[, 1], nvx), nrow = nsx, byrow = F)
    A1=(buffer - matrix(colMeans(buffer), nrow = nsx,ncol = ncol(buffer), byrow = T))
    B=(X - matrix(colMeans(X), nrow = nsx, ncol = ncol(X), byrow = TRUE))
    D=A1*B
    G=colMeans(A1*B)
    CC=G/(nsx - 1) * nsx
    s1 <- sd(model$T[, 1])
    S1 <- matrix(rep(s1, nvx),nrow = 1, byrow = T)
    S2 <- apply(X, 2, sd)
    median_R <- CC / (S1 * S2)
    median_P <- tTest(median_R, nsx)# pearson correlation significant test
    if (!is.null(pCutoff)) {
      IX_P <- which(median_P < pCutoff)#
    }


    ## save out put
    model$sig_idx <- IX_P  # the significant variables
    model$signs <- Signs   # up- or down-regulate

  }
  else {
    # CV
    # nrcv-fold cross validation
    block_num <- floor(nsx / nrcv)
    Q2Yhatcum <- rep(0, nrcv)


    for (cv in 1:nrcv){# nrcv iterations of CV
      idx_test <- ((1:block_num) - 1) * nrcv + cv
      idx_tr <- 1:nsx
      idx_tr <- idx_tr[!(idx_tr %in% idx_test)]
      X_test <- X[idx_test, ]
      Y_test <- Y[idx_test, ]
      X_tr <- X[idx_tr, ]
      Y_tr <- Y[idx_tr, ]



      nsx_tr <- nrow(X_tr)
      nvx_tr <- ncol(X_tr)
      nsy_tr <- nrow(Y_tr)
      nvy_tr <- ncol(Y_tr)
      nsx_test <- nrow(X_test)
      nvx_test <- ncol(X_test)
      nsy_test <- nrow(Y_test)
      nvy_test <- ncol(Y_test)


      #Univariance scaling
      X_tr <- scale(X_tr)
      Y_tr <- scale(Y_tr)
      X_test <- scale(X_test)
      Y_test <- scale(Y_test)
      X_test <- scale(X_test, center = colMeans(X_test, na.rm = TRUE), scale = apply(X_test, 2, sd, na.rm = TRUE))
      Y_test <- scale(Y_test, center = colMeans(Y_test, na.rm = TRUE), scale = apply(Y_test, 2, sd, na.rm = TRUE))

      #training
   
      model <- mjrO2pls(X_tr, Y_tr, nc, ncox, ncoy, "standard")
      
       
      # train OPLS model
      #yhat
      modelPredy=mjrO2plsPred(X_test, Y_test, model, ncox, ncoy, 'x')# prediction using OPLS model
      # #the overall Q2
      SSY <- sum(Y_test^2)
      Q2Yhatcum[cv] <- 1 - sum((modelPredy$Yhat - Y_test)^2) / SSY

    }
    model$Q2Yhatcum <- mean(Q2Yhatcum[!is.na(Q2Yhatcum)], na.rm = TRUE)
  }
   
  return(model)

}

