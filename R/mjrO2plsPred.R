#' Used to predict the results of the test set data X_test and Y_test in the X direction, use the trained O2PLS model model and specify the number of orthogonal components ncox in the X direction and the number of orthogonal components ncoy in the Y direction.
#'
#' @param X X data for test set.
#' @param Y Y data for test set.
#' @param model Stores the relevant parameters and weights of the O2PLS model.
#' @param oax Number of orthogonal components in the X direction.
#' @param oay Number of orthogonal components in the Y direction.
#' @param dir dir: direction of prediction, can be 'x' or 'y'.
#'
#' @return modelPred,The score matrix T, Y prediction matrix Yhat and orthogonal prediction To in the X direction are stored in the modelPred structure and returned.


mjrO2plsPred <- function(X, Y, model, oax, oay, dir) {
  #------------------------------------------------
  #Author: Mattias Rantalainen, Imperial College, 2007
  #Copyright: Mattias Rantalainen, 2007
  #------------------------------------------------
  if(tolower(dir) == 'x') {
    To <- list()
    #if(oax > 0)

    for(i in 1:oax) {
      to <- X %*% model$Wo[, i] %*% solve(t(model$Wo[, i]) %*% model$Wo[, i])
      To[[i]] <- to
      X <- X - to %*% t(model$Pyo[, i])
    }


    T <- X %*% model$W %*% solve(t(model$W) %*% model$W)
    Yhat <- T %*% model$Bts[[oax + 1, oay + 1]] %*% t(model$C)
    modelPred <- list(T = T, Yhat = Yhat, To = To)

    #}
  }

  if(tolower(dir) == 'y') {
    Uo <- list()
    #  if(oay > 0)
    for(i in 1:oay) {
      uo <- Y %*% model$Co[, i] %*% solve(t(model$Co[, i]) %*% model$Co[, i])
      Uo <- cbind(Uo, uo)
      Y <- Y - uo %*% t(model$Pxo[, i])
    }


    U <- Y %*% model$C %*% solve(t(model$C) %*% model$C)
    Xhat <- U %*% model$Bus[[oax + 1, oay + 1]] %*% t(model$W)
    modelPred <- list(U = U, Xhat = Xhat, Uo = Uo)

    #}
  }
  return(modelPred)
}

