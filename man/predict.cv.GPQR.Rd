\name{predict.cv.GPQR}
\alias{predict.cv.GPQR}
\title{make predictions from a "cv.GPQR" object.}
\description{
This function makes predictions from a cross-validated \code{GPQR} model,
using the stored \code{"GPQR"} object, and the optimal value
chosen for \code{lambda}.
}

\usage{
predict.cv.GPQR(cv, newx, s = c("lambda.1se", "lambda.min"), ...)
}


\arguments{

\item{cv}{fitted \code{\link{cv.GPQR}} object.}
\item{newx}{matrix of new values for \code{x} at which predictions are to be made. Must be a matrix. See documentation for \code{predict.GPQR}.}
\item{s}{value(s) of the penalty parameter \code{lambda} at which predictions are required. Default is the value \code{s="lambda.1se"} stored on the CV object. Alternatively \code{s="lambda.min"} can be used}
\item{\dots}{not used. Other arguments to predict. } 

}

\author{Mohamed Ouhourane, Yi Yang, Karim Oualkacha
Maintainer: Mohamed OUHOURANE <mohamed.ouhourane@gmail.com>}
\references{
 Mohamed Ouhourane, Yi Yang, Karim Oulkacha: "Group Penalized Smooth Quantile Regression" 
}

\examples{
data(bardet)
group <- rep(1:20,each=5)
cv <- cv.GPQR(x=bardet$x,y=bardet$y,group=group,method="GLasso",check="f1",taux=0.5)
predict.cv.GPQR(cv,bardet$x,s = "lambda.1se")
}

