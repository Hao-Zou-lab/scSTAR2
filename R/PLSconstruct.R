
#' function PLSconstruct constructs PLS model based on the given data:data1 and data2, and parameter specifications.
#'
#' @param data1 the data1 matrix from control and case groups with genes by cells.
#' @param data2 the data2 matrix from control and case groups with genes by cells.
#' @param prep preprocessing methods.
#' @param NCV number of folds of cross-validation.
#' @param PLScomp_def the number of PLS components, 0 - the value is automatically estimated;otherwise used the given value.
#' @param minNC the minimum number of PLS components.
#'
#' @return MODEL,the PLS model containing all PLS parameters and outputs



## PLS1 projection
PLSconstruct <- function(data1, data2,prep, NCV, PLScomp_def, minNC)
{
  # function PLSconstruct constructs PLS model based on the given data: data1 and data2,
  #         and parameter specifications.
  #
  # Input: data1 and data2 - the data matrix from control and case groups with genes by cells.
  #        prep - preprocessing methods.
  #        NCV  - number of folds of cross-validation.
  #        PLScomp_def - the number of PLS components, 0 - the value is automatically estimated;
  #                   otherwise used the given value.
  #        minNC - the minimum number of PLS components.
  # Output: model - the PLS model containing all PLS parameters and outputs

  # reshape data for PLSDA model construction
  X = rbind(data1, data2)
  Y = matrix(0, dim(X)[1],2)
  Y[seq_len(dim(data1)[1]),1] <- 1
  Y[(dim(data1)[1]+1):dim(Y)[1],2] <- 1

  maxNC = 5

  # scaling the data
  tmp <-dim(X)
  nsx <- tmp[1]
  nvx <- tmp[2]
  tmp <- dim(Y)
  nsy <- tmp[1]
  nvy <- tmp[2]

  muX = colMeans(X)
  muY = colMeans(Y)

  sigmaX = apply(X,2,sd)
  sigmaY = apply(Y,2,sd)

  model <- list(NULL)

  if (prep == 'no')
  {
    model$preprocessing<-'no'
  } else if (prep == 'mc')
  {
    #disp('Meancentering')
    model$preprocessing<-'mc'
    X<-X-matrix(rep(muX,each = nsx),nsx)
    Y<-Y-matrix(rep(muY,each = nsy),nsy)
  } else if (prep == 'uv')
  {
    # disp('Univariance scaling')
    model$preprocessing<-'uv'
    X<-X-matrix(rep(muX,each = nsx),nsx)
    Y<-Y-matrix(rep(muY,each = nsy),nsy)
    X<-X/matrix(rep(sigmaX,each = nsx),nsx)
    Y<-Y/matrix(rep(sigmaY,each = nsy),nsy)
  } else if (prep == 'pa')
  {
    model$preprocessing<-'pa'
    X<-X-matrix(rep(muX,each = nsx),nsx)
    Y<-Y-matrix(rep(muY,each = nsy),nsy)
    X<-X/matrix(rep(sqrt(sigmaX),each = nsx),nsx)
    Y<-Y/matrix(rep(sqrt(sigmaY),each = nsy),nsy)
  } else
  {
    model<-NULL
    stop('Unknown Preprocessing\n')
  }


  # PLS model construction
  block_num = floor(dim(X)[1]/NCV)
  Q2_ori_1 = 0
  Q2_ori_2 = 0
  nsx = dim(X)[1]
  if (PLScomp_def == 0)  # estimate the optimal number of PLScomp
  {
    PLScomp = 1
    for (nc in 1:maxNC)
    {
      Q2Yhatcum = matrix(0,1,NCV)
      for (cv in 1:NCV)
      {
        # cross validation
        idx_test = NCV*(1:block_num-1)+cv
        idx_tr = 1:nsx
        idx_tr = idx_tr[-idx_test]
        X_test = X[idx_test,]
        Y_test = Y[idx_test,]
        X_tr = X[idx_tr,]
        Y_tr = Y[idx_tr,]

        # model construct
        plsrDT = plsr(Y_tr ~ X_tr,nc)

        BETA <- plsrDT$coefficients[,,nc, drop=FALSE]
        # cumulative coefficients Intercept = T
        dB <- dim(BETA)
        dB[1] <- dB[1] + 1
        dnB <- dimnames(BETA)
        dnB[[1]] <- c("(Intercept)", dnB[[1]])
        BInt <- array(dim = dB, dimnames = dnB)
        BInt[-1,,] <-  BETA

        for (i in seq(along = nc))
          BInt[1,,i] <- plsrDT$Ymeans - plsrDT$Xmeans %*% BETA[,,i]

        BETA <- BInt[,,1]

        # predicting on the model
        PredY = cbind(matrix(1,dim(X_test)[1],1),X_test)%*%BETA
        # calculating Q2
        SSY = sum(Y_test^2)
        Q2Yhatcum[cv] = 1-sum((PredY-Y_test)*(PredY-Y_test))/SSY
      }
      Q2 = mean(Q2Yhatcum[!is.na(Q2Yhatcum)])

      # check if the current Q2 value is maximum
      if (Q2_ori_2>Q2 & Q2_ori_2>Q2_ori_1 & Q2_ori_1>Q2)
      {
        PLScomp = nc - 2  # this is the critical number of PLS components
        break
      } else {
        Q2_ori_2 = Q2_ori_1
        Q2_ori_1 = Q2
      }
    }
    PLScomp = max(PLScomp, minNC)
  } else
  {
    PLScomp = PLScomp_def
  }

  # model construction with the estimated number of PLS components
  plsestcom = plsr(Y ~ X ,PLScomp)
  XL = loadings(plsestcom)

  BETA <- plsestcom$coefficients[,,PLScomp, drop=FALSE]
  # cumulative coefficients Intercept = T
  dB <- dim(BETA)
  dB[1] <- dB[1] + 1
  dnB <- dimnames(BETA)
  dnB[[1]] <- c("(Intercept)", dnB[[1]])
  BInt <- array(dim = dB, dimnames = dnB)
  BInt[-1,,] <-  BETA

  for (i in seq(along = PLScomp))
    BInt[1,,i] <- plsestcom$Ymeans - plsestcom$Xmeans %*% BETA[,,i]

  BETA <- BInt[,,1]


  model$XL = matrix(XL, dim(XL)[1],dim(XL)[2])
  model$BETA = BETA

  if (prep == 'no')
  {
  } else if (prep == 'mc')
  {
    model$mu = muX
  } else if (prep == 'uv')
  {
    model$mu = muX
    model$sigmaX = sigmaX
  } else if (prep == 'pa')
  {
    model$mu = muX
    model$sigmaX = sigmaX
  }

  return(model)
}
