


#' Plot correlations between baselearners in a prediction rule ensemble (pre)
#' 
#' \code{corplot} plots correlations between baselearners in a prediction rule ensemble
#'  
#' @param object object of class pre
#' @inheritParams print.pre
#' @param colors vector of contiguous colors to be used for plotting. If 
#' \code{colors = NULL} (default), \code{colorRampPalette} is used to generate
#' a sequence of 200 colors going from red to white to blue. A different set of 
#' plotting colors can be specified here, for example: 
#' \code{cm.colors(100)}, \code{rainbow_hcl)(100)} (the latter requires
#' package \code{colorspace}). 
#' or \code{colorRampPalette(c("red", "yellow", "green"))(100)}.
#' @param fig.plot plotting region to be used for correlation plot. See 
#' \code{fig} under \code{\link{par}}.
#' @param fig.legend plotting region to be used for legend. See \code{fig} 
#' under \code{\link{par}}.
#' @param legend.breaks numeric vector of breakpoints to be depicted in the 
#' plot's legend. Should be a sequence from -1 to 1.
#' @examples \donttest{set.seed(42)
#' airq.ens <- pre(Ozone ~ ., data = airquality[complete.cases(airquality),])
#' corplot(airq.ens)}
#' @export
corplot <- function(object, penalty.par.val = "lambda.1se", colors = NULL,
                    fig.plot = c(0, 0.85, 0, 1), fig.legend = c(.8, .95, 0, 1),
                    legend.breaks = seq(-1, 1, by = .1)) {
  
  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(
      c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", 
        "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(200)
  }
  ## get coefficients:
  coefs <- as(coef(object$glmnet.fit, s = penalty.par.val), Class = "matrix")
  ## create correlation matrix of non-zero coefficient learners:
  cormat <- cor(object$modmat[,names(coefs[coefs != 0,])[-1]])
  ## set layout for plotting correlation matrix and legend::
  layout(matrix(rep(c(1, 1, 1, 1, 2), times = 5), 5, 5, byrow = TRUE))
  ## set region for correlation matrix:
  par(fig = fig.plot)
  ## plot correlation matrix:
  image(x = 1:nrow(cormat), y = 1:ncol(cormat), z = unlist(cormat), 
        axes = FALSE, zlim = c(-1, 1),  
        xlab = "", ylab = "", srt = 45, 
        col = colors)
  axis(1, at = 1:nrow(cormat), labels = colnames(cormat), las = 2)
  axis(2, at = 1:ncol(cormat), labels = colnames(cormat), las = 2)
  ## plot legend:
  par(fig = fig.legend, new = TRUE)
  col_inds <- round(seq(1, length(colors), length.out = length(legend.breaks)-1))
  image_scale(z = legend.breaks,
              col = colors[col_inds],
              axis.pos = 4,
              breaks = legend.breaks)
  axis(4, at = legend.breaks, las = 2)
  return(invisible(cormat))
}


## Internal function for creating legend for corplot:
image_scale <- function(z, col, breaks, axis.pos = 4, add.axis = TRUE) {
  poly <- vector(mode = "list", length(col))
  for(i in seq(poly)) {
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if (axis.pos %in% c(1,3)) { 
    ylim <- c(0,1)
    xlim <- range(breaks)
  }
  if (axis.pos %in% c(2,4)){
    ylim <- range(breaks)
    xlim <- c(0,1)
  }
  plot(1, 1, t = "n", ylim = ylim, xlim = xlim, axes = FALSE, xlab = "", 
       ylab = "", xaxs = "i", yaxs = "i")  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}
