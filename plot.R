#' Plot the estimated eigenfunctions
#'
#' Plot the estimated eigenfunctions from \code{LVPCA} outputs.
#'
#' @name plot.LVPCA
#' @aliases plot.LVPCA
#' @param x An output from \code{LVPCA}
#' @param k The number of principal components to plot.
#' @param type The type of plot with default \code{"b"}.
#' @param xlab A title for the x axis with default \code{"m"}.
#' @param ylab A title for the y axis with default \code{""}.
#' @param cex.lab The cex for lab with default \code{3}.
#' @param cex.axis The cex for axis with default \code{2.5}.
#' @param cex.main The cex for main with default \code{3.5}
#' @param xaxt A character specifies the x axis type with default \code{"n"}.
#'
#' @import scales
#' @importFrom graphics abline axis par plot title
#' @export
#' @details If \code{k} is not specified, then \code{k} from the \code{LVPCA} output is used. If \code{model = "2WayNested"} and
#' \code{k} is different between two levels, the smaller number between \code{k.w} and \code{k.z} is used.
#' @return
#' If argument \code{model = "1Way"}, then the top \code{k} eigenfunctions are plotted, with titles showing the fraction of variance explained
#' by each component.
#' If argument \code{model="2WayNested"}, the top \code{k} subject level eigenfunctions are plotted in the first row and the top \code{k}
#' replicated-within-subject level eigenfunctions are plotted in the second row, with titles showing the fraction of variance explained out of
#' the level-specific total variance.
#' @examples
#' data("data.2levels")
#' result1=LVPCA(data.2levels,m=3,L1=100,L2=5, FVE_threshold = 0.9, model="2WayNested")
#' plot(result1)

#require(scales)
#library(scales)

plot.LVPCA = function(x, k=NULL, type="b", cex.main=3.5,
                      ylab = "", cex.axis = 2.5, xlab = "m", cex.lab=3, xaxt = 'n'){

  #if (missing(y)) {y=x}
  result=x
  m=result$option$m
  p=ncol(result$predx)/m
  if (result$option$model=="1Way"){
    if (is.null(k)){
       k=result$option$k
    }
    par(mfrow=c(1,k), oma=c(0.5,0.5,0.5,0.5), mar = c(5.5,4.5,5,2), mgp=c(3.5,1.5,0))
    for (i in 1:k){
      plot(result$phi[,i], type=type, ylab = ylab, cex.axis = cex.axis, xlab = xlab, cex.lab=cex.lab, xaxt = xaxt)
      abline(v=seq(0,m*p,by=p), col=rep("black", 4),lwd=1)
      axis(side=1,at=(0.5+0:(m-1))*p,labels=1:m, cex.axis=cex.axis)
      title(bquote(hat(phi)[.(i)]~"("~.(percent(result$FVE[i]))~")"),cex.main=cex.main,line=2.3)
    }


  }else if (result$option$model=="2WayNested"){
    if (is.null(k)){
       k=min(result$option$k$k.z,result$option$k$k.w)
    }
    par(mfrow=c(2,k), oma=c(0.5,0.5,0.5,0.5), mar = c(5.5,4.5,5,2), mgp=c(3.5,1.5,0))
    for (i in 1:k){
      plot(result$phi.z[,i], type=type, ylab = ylab, cex.axis = cex.axis, xlab = xlab, cex.lab=cex.lab, xaxt = xaxt)
      abline(v=seq(0,m*p,by=p), col=rep("black", 4),lwd=1)
      axis(side=1,at=(0.5+0:(m-1))*p,labels=1:m, cex.axis=cex.axis)
      title(bquote(hat(phi)[.(i)]^"z"~"("~.(percent(result$FVE.z[i]))~")"),cex.main=cex.main,line=2.3)
    }
    for (i in 1:k){
      plot(result$phi.w[,i], type=type, ylab = ylab, cex.axis = cex.axis, xlab = xlab, cex.lab=cex.lab, xaxt = xaxt)
      abline(v=seq(0,m*p,by=p), col=rep("black", 4),lwd=1)
      axis(side=1,at=(0.5+0:(m-1))*p,labels=1:m, cex.axis=cex.axis)
      title(bquote(hat(phi)[.(i)]^"w"~"("~.(percent(result$FVE.w[i]))~")"),cex.main=cex.main,line=2.3)
    }

  }


}
