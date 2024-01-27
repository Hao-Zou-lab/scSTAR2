#' This function implements model building and parameter estimation for the O2PLS method and provides a list containing model information. Such models can be used to further predict and explain multivariate data.
#'
#' @param X Data for explanatory variables.
#' @param Y Data for the response variable.
#' @param pax Number of principal components in the X direction.
#' @param oax Number of orthogonal components in the X direction.
#' @param oay Number of orthogonal components in the Y direction.
#' @param mtype Select the calculation method, which can be 'standard' or 'large'.
#' @param varargin Optional parameter, mainly used to specify the block size when calculating large matrices.
#'
#' @return model,A structure containing the trained O2PLS model and related statistical information.


mjrO2pls = function(X,Y,pax,oax,oay,mtype,varargin){
  #------------------------------------------------
  #Author: Mattias Rantalainen, Imperial College, 2007
  #Copyright: Mattias Rantalainen, 2007
  #------------------------------------------------
  varargin = NULL
  if (!is.null(varargin) && length(varargin) > 0 && !is.null(varargin[[1]])) {#double checks here for backward compatibilty
    splitSize <- varargin[[1]]
  } else {
    splitSize <- 300  # this is the Y block size (if mtype='large')
  }

  ssx = sum(colSums(X^2))
  ssy = sum(colSums(Y^2))
  #large matrices
  if(mtype == "large"){
    m = nrow(Y)
    n = ncol(Y)
    nsplits = floor(n/splitSize)
    Wtot = matrix()
    Ctot = matrix()

    h = tcltk::tkProgressBar(title = "Progress", label = "estimating Wp for Y-blocks...", min = 0,max = 100)#,num2str(1),' of ',num2str(nsplits)])
    for (i in nsplits:1) {
      #disp(['ysplit',num2str(i)])
      tcltk::setTkProgressBar(h, value = (i-1)/nsplits)#,['estimating Wp for Y-block: ',num2str(i),' of ',num2str(nsplits)])
      svd_result = svd(t(Y[, ((i-1)*splitSize+1):(i*splitSize)]) %*% X)
      Wtmp = svd_result$u
      Stmp = diag(svd_result$d)
      Ctmp = svd_result$v
      Wtot <- cbind(Wtot, Wtmp %*% Stmp)


      if (i == nsplits && (n %% splitSize) != 0) {
        #disp(['last Y-block - got remaining part of ',num2str(mod(n,splitSize)),' vars...']);
        svd_resultt <- svd(t(Y[, ((i * splitSize) + 1):ncol(Y)]) %*% t(X))
        Wtmp <- svd_resultt$u
        Stmp <- diag(svd_resultt$d)
        Ctmp <- svd_resultt$v
        if (mod(n %% splitSize, splitSize) < pax) {
          paxTmp <- mod(n %% splitSize, splitSize)
        }
        else {
          paxTmp <- pax
        }

        Wtmp <- svd_resultt$u[, 1:paxTmp]
        Stmp <- diag(svd_resultt$d[1:paxTmp])
        Wtot <- cbind(Wtot, Wtmp %*% Stmp)
      }
    }

    tcltk::tkdestroy(h)


    svd_result_1 <- svd(Wtot)
    Wtot <- NULL
    W <- svd_result_1$u[, 1:pax]
    S <- diag(svd_result_1$d[1:pax])
    Ts [[1]]<- list(X %*% W)
    T <- Ts[[1]]
    C <- t(Y) %*% T
    C <- C %*% diag(1 / sqrt(colSums(C[, 1:pax]^2)))
  }
  else { #conventional
    svd_result_2 <- svd(t(Y) %*% X)
    C <- svd_result_2$u
    W <- svd_result_2$v
    S <- diag(svd_result_2$d)
    C = C[, 1:pax]
    W = W[, 1:pax]
    S = S[1:pax, 1:pax]

    #     [~,~,~,~,W,~,C,~] = plsr(X,Y,pax)
    Ts <- list(X %*% W)#
    T <- Ts[[1]]
  }

  Exy <- X - T %*% t(W)
  Xr <- X

  R2Xo <- 0
  R2Xcorr <- sum(sum((Ts[[1]] %*% t(W))^2)) / ssx
  R2X <- 1 - sum(sum(Exy^2)) / ssx


  Wo <- vector()
  Pyo <- vector()
  To <- vector()

  for (i in 1:oax) {
    svd_result_3 <- svd(t(Exy) %*% Ts[[i]])
    wo <- svd_result_3$u[, 1] #/sqrt(wo'*wo)
    syo <- diag(svd_result_3$d)
    wor <- svd_result_3$v
    to <- Xr %*% wo  #/(wo'*wo)
    pyo <- t(Xr) %*% to / as.numeric(t(to) %*% to)
    Xr <- Xr - to %*% t(pyo)

    Wo <- cbind(Wo, wo)
    Pyo <- cbind(Pyo, pyo)
    To <- cbind(To, to)

    Ts[[i + 1]] <- Xr %*% W
    Exy <- X - Ts[[i + 1]] %*% t(W) - To %*% t(Pyo)

    R2Xo <- c(R2Xo, sum(sum((To %*% t(Pyo))^2)) / ssx)
    R2Xcorr <- c(R2Xcorr, sum(sum((Ts[[i + 1]] %*% t(W))^2)) / ssx)
    R2X <- c(R2X, 1 - sum(sum(Exy^2)) / ssx)
    T <- Ts[[i + 1]]
  }



  Yr <- Y
  Us <- list()
  Us[[1]] <- Yr %*% C
  U <- Us[[1]]
  Fxy <- Y - Us[[1]] %*% t(C)


  R2Yo <- 0
  R2Ycorr <- sum(sum((Us[[1]] %*% t(C))^2)) / ssy
  R2Y <- 1 - sum(sum(Fxy^2)) / ssy

  Uo <-vector()
  Pxo <- vector()
  Co <- vector()

  if (oay > 0) {
    for (i in 1:oay) {
      svd_result_4 <- svd(t(Fxy) %*% Us[[i]], nu = 0)
      co <- svd_result_4$u[, 1] #/sqrt(co'*co)
      sxo <- svd_result_4$d #/(co'*co);
      cor <- svd_result_4$v
      uo <- Yr %*% co
      pxo <- t(Yr) %*% uo /(t(uo) %*% uo)
      Yr <- Yr - uo %*% t(pxo)

      Co <- cbind(Co, co)
      Pxo <- cbind(Pxo, pxo)
      Uo <- cbind(Uo, uo)

      Us[[i + 1]] <- Yr %*% C
      Fxy <- Y - Us[[i + 1]] %*% t(C) - Uo %*% t(Pxo)

      R2Yo <- c(R2Yo, sum(sum((Uo %*% t(Pxo))^2)) / ssy)
      R2Ycorr <- c(R2Ycorr, sum(sum((Us[[i + 1]] %*% t(C))^2)) / ssy)
      R2Y <- c(R2Y, 1 - sum(sum(Fxy^2)) / ssy)
      U <- Us[[i + 1]]
    }
  }
  Bus <- matrix(list(), nrow = oax + 1, ncol = oay + 1)
  Bts <- matrix(list(), nrow = oax + 1, ncol = oay + 1)
  for (i in 1:(oax + 1)) {
    for (j in 1:(oay + 1)) {
      Bus[[i, j]] <- solve(t(Us[[j]]) %*% Us[[j]]) %*% t(Us[[j]]) %*% Ts[[i]]
      Bts[[i, j]] <- solve(t(Ts[[i]]) %*% Ts[[i]]) %*% t(Ts[[i]]) %*% Us[[j]]
    }
  }

  Bu <- solve(t(U) %*% U) %*% t(U) %*% T
  Bt <- solve(t(T) %*% T) %*% t(T) %*% U

  R2Yhat <- vector()
  R2Xhat <- vector()



  #This is OC style - keeping for now... :
  for (i in 1:(oax + 1)) {
    BtTmp <- solve(t(Ts[[i]]) %*% Ts[[i]]) %*% t(Ts[[i]]) %*% U
    YhatTmp <- Ts[[i]] %*% BtTmp %*% t(C)
    R2Yhat <- c(R2Yhat, 1 - sum(sum((YhatTmp - Y)^2)) / ssy)
  }
  for (i in 1:(oay + 1)) {

    BuTmp <- solve(t(Us[[i]]) %*% Us[[i]]) %*% t(Us[[i]]) %*% T
    XhatTmp <- Us[[i]] %*% BuTmp %*% t(W)
    R2Xhat <- c(R2Xhat, 1 - sum(sum((XhatTmp - X)^2)) / ssx)
  }
  model <- list(
    T = T,
    Ts = Ts,
    W = W,
    Wo = Wo,
    Pyo = Pyo,
    To = To,
    U = U,
    Us = Us,
    C = C,
    Co = Co,
    Pxo = Pxo,
    Uo = Uo,
    Bt = Bt,
    Bu = Bu,
    Bts = Bts,
    Bus = Bus,
    R2X = R2X,
    R2Xcorr = R2Xcorr,
    R2Xo = R2Xo,
    R2Xhat = R2Xhat,
    R2Y = R2Y,
    R2Ycorr = R2Ycorr,
    R2Yo = R2Yo,
    R2Yhat = R2Yhat,
    ssx = ssx,
    ssy = ssy
  )
}
