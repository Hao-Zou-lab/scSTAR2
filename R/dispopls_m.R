#'Display the result of an opls modelling.
#'The first plot is the first and second PLS scores of the two classes.
#'The second plot is the weights of latent variables.
#'
#' @param m opls model (see O2PLS).
#' @param X Explicative variable .
#' @param Y Response variables (dummy variables in case of Discriminant analysis).
#' @param xscale X matrix unit scale (to plot the opls coefficients with the proper axis e.g. cellindex).
#'
#' @return Sign_buffer,For some elements in the buffer, if they are statistically significant and they are upregulated relative to the baseline, then the corresponding element values in sign_buffer will be positive; if they are significant and relative is lower than the baseline, then the element value in the corresponding sign_buffer will be negative.
#' @return p2, If there are significant variables with different positive and negative signs, mark them with red and green scatter points.
#' @return p3, If there are only significant variables with positive sign, mark them with red scatter points.
#' @return p4, If there are only significant variables with negative sign, mark them with green scatter points.
#' @return p5, If there is no significant variable, keep the original histogram.


dispopls_m <- function(m, X, Y, xscale){
  #
  # Example of command: dispopls_m(model,X,Y,ppm, {'1','2'}, 'oplsda for X vs Y', 1, 1);
  #
  # Display the result of an opls modelling.
  # The first plot is the first and second PLS scores of the two classes.
  # The second plot is the weights of latent variables.
  #
  # Input:
  # m      - opls model (see O2PLS),
  # X      - Explicative variable (ex: NMR spectra, MS data),
  # Y      - Response variables (dummy variables in case of Discriminant analysis)
  # xscale - X matrix unit scale (to plot the opls coefficients
  #                 with the proper axis e.g.cellindex)
  # DIET_list - class name e.g. for class 1 vs class 2 put {'1','2'}
  # figTitle - title of the figure e.g. for 'oplsda for X vs Y'
  # flag   - '1' plot spectra, '0' not
  # plotOp - '1' plot correlated weights, '0' orthogonal weights
  # cutoff - minimum correlation to be displayed


  #
  # Output:
  # Correlations - coefficient de correlation between the response
  #                       variables and the explicative variables
  # Covariances  - covariances between the response variables
  #                       and the explicative variables
  # Xin Zou,
  height <- 0.7
  buffer <- apply(X[Y[, 2] == 1, ], 2, median) - apply(X[Y[, 1] == 1, ], 2, median)
  buffer_1<-max(buffer)
  if(buffer_1<=0){buffer_1=min(buffer)}
  xscale <- c(1:ncol(X))
  data_to_plot <- data.frame(xscale = xscale, buffer = buffer)
  #x11()
  p1 <- ggplot(data_to_plot, aes(x = xscale, y = buffer)) +
    geom_bar( stat = "identity",fill= "blue", color = "black") +
    labs(x = "cellindex", y = "O-PLS coefficients (a.u.)") +
    theme( plot.margin = margin(0.12, 1.0-((height+0.12)) ,0.75,height),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = "black"),
           axis.line.x=element_line(linetype=1,color="black",size=1),
           axis.line.y=element_line(linetype=1,color="black",size=1),
           axis.ticks.y=element_line(color="black",size=1,lineend = 10),
           axis.ticks.x=element_line(color="black",size=1,lineend = 1),
           axis.text = element_text(size = 14),
           axis.title = element_text(size = 14),
           panel.border = element_rect(color = "black", fill = NA, size = 1.2))
  print(p1)

  ## significant variables
  idx <- m$sig_idx
  buffer <- apply(X[Y[, 2] == 1, ], 2, median) - apply(X[Y[, 1] == 1, ], 2, median)

  sign_buffer <- sign(buffer[idx])


  ## colour coded accoding to regulation trend, up-regulated compared to
  # baseline 'red'; down-regulate 'orange'
  if (length(which(sign_buffer < 0)) > 0 && length(which(sign_buffer > 0)) > 0) {
    positive_data <- data.frame(x = xscale[idx[which(sign_buffer > 0)]],y=0)
    negative_data <- data.frame(x=xscale[idx[which(sign_buffer < 0)]],y=0)
    #x11()
    p2 <-p1 +
      geom_point(data = positive_data, aes(x = x, y = y), color = 'red', size = 2) +
      geom_point(data = negative_data, aes(x = x, y = y), color = 'green', size = 2) +
      theme(plot.margin = margin(0.12, 0.88, 0.75, 0.020))
    print(p2)

  }else if (length(which(sign_buffer > 0)) > 0 && length(which(sign_buffer < 0)) == 0){
    positive_data <- data.frame(x = xscale[idx[which(sign_buffer > 0)]],y=0)
    #x11()
    p3 <-p1 +
      geom_point(data = positive_data, aes(x = x, y = y), color = 'red', size = 2) +
      theme(plot.margin = margin(0.12, 0.88, 0.75, 0.020))
    print(p3)

  }else if(length(which(sign_buffer < 0)) > 0 && length(which(sign_buffer > 0)) == 0){
    negative_data <- data.frame(x=xscale[idx[which(sign_buffer < 0)]],y=0)
    #x11()
    p4 <-p1 +
      geom_point(data = negative_data, aes(x = x, y = y), color = 'green', size = 2) +
      theme(plot.margin = margin(0.12, 0.88, 0.75, 0.020))
    print(p4)

  }else{
    #x11()
    p5 <-p1
    print(p5)

  }
  model_1 <- list()
  model_1$Sign_buffer=sign_buffer
  model_1$buffer_1=buffer_1
  if (length(which(sign_buffer < 0)) > 0 && length(which(sign_buffer > 0)) > 0){
    model_1$p2 = p2
  } else if (length(which(sign_buffer > 0)) > 0 && length(which(sign_buffer < 0)) == 0){
    model_1$p3 = p3
  } else if (length(which(sign_buffer < 0)) > 0 && length(which(sign_buffer > 0)) == 0){
    model_1$p4 = p4
  } else {
    model_1$p5 = p5
  }

  return(model_1)

}
