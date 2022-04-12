# Get lambda.min and lambda.1se
getmin <-function(lambda, cvm, cvsd) {
cvmin <- min(cvm)
idmin <- cvm <- cvmin
lambda.min <- max(lambda[idmin])
idmin<- match(lambda.min, lambda)
semin<-(cvm + cvsd)[idmin]
idmin <- cvm <- semin
lambda.1se <- max(lambda[idmin])
list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}


#This function produce a coefficient profile plot of the coefficient
#paths for a fitted GPQR object.

GPQR_illustration <-function(penalty, taux){
  cv <- cv.GPQR(x=Xtr,y=Ytr,group=group,method=penalty,check="f1",taux=taux)
  fittGL=t(cv$finalfit$beta)
  seqLambda=cv$lambda
  seqLambda=cv$lambda
  ll = getmin(cv$lambda, cv$cv, cv$cv.error)
  l = which(cv$lambda == ll$lambda.min)
  main_lab = paste("Q-",penalty," ",expression(tau)," = ", taux)
  matplot(seqLambda,fittGL[,cc1], type = "l",col = 3,lty = 1,ylim=xlm,lwd=1,
          ylab="Coefficients",main=main_lab, xlab = expression(lambda))
  matlines(seqLambda,fittGL[,cc2], type = "l",col = 2,lty = 1,lwd=1)
  matlines(seqLambda,fittGL[,cc3], type = "l",col = 4,lty = 1,lwd=1)
  matlines(seqLambda,fittGL[,cc4], type = "l",col = 1,lty = 1,lwd=1)
  matlines(seqLambda,fittGL[,20], type = "l",col = 6,lty = 1,lwd=1)
  abline(v=seqLambda[l],col=1,lty = 2)
  text(0.02,0.75,expression("G"[11]),col=6)
  text(0.02,3.2,expression("G"[1]),col=3)
  text(0.02,2.2,expression("G"[2]),col=2)
  text(0.02,-1.4,expression("G"[3]),col=4)
}
# This function produce a coefficient profile plot of the coefficient paths for
#a fitted grpreg object.
grpreg_illustration <-function(penalty){
  cv <- cv.grpreg(Xtr, Ytr, group, penalty=penalty)
  fittL=t(cv$fit$beta)[,-1]
  seqLambdaL=cv$lambda
  l=cv$min
  matplot(seqLambdaL,fittL[,cc1], type = "l",col = 3,lty = 1,ylim=xlm,lwd=1,
          ylab="Coefficients",main=paste("LS-",penalty) ,xlab = expression(lambda))
  matlines(seqLambdaL,fittL[,cc2], type = "l",col = 2,lty = 1,lwd=1)
  matlines(seqLambdaL,fittL[,cc3], type = "l",col = 4,lty = 1,lwd=1)
  matlines(seqLambdaL,fittL[,cc4], type = "l",col = 1,lty = 1,lwd=1)
  matlines(seqLambdaL,fittL[,20], type = "l",col = 6,lty = 1,lwd=1)
  ###abline(v=seqLambdaL[dd],col = 1,lty = 3)
  ##abline(v=seqLambdaL[dd+ddd],col = 2,lty = 3)
  abline(v=seqLambdaL[l],col=1,lty = 2)
  text(0.4,0.75,expression("G"[11]),col=6)
  text(0.4,3.3,expression("G"[1]),col=3)
  text(0.4,2.3,expression("G"[2]),col=2)
  text(0.4,-1.5,expression("G"[3]),col=4)
}
# The function "plot_grpreg" produce the optimal value  for the regression
#coefficients of the grpreg-methods are shown as a function of the genomic position.
plot_grpreg <- function (penalty){
  cvfit = cv.grpreg(Xc, pheno, group, penalty=penalty)
  coefGM=predict(cvfit, Xc, type="coefficients")[-1]
  matplot(positionC/1e6,coefGM,col=4, type = "p",pch=19,cex=1,ylab=ylab,xlab="Location",main=paste("LS-",penalty),ylim=c(-0.008,0.006))
  abline(v=11.235,lty=2)
  abline(v=11.385,lty=2)
}
#The function "plot_GPQR" produce the optimal value  for the regression coefficients of
#the GPQR-methods are shown as a function of the genomic position.
plot_GPQR <- function (penalty,taux){
  cvGS=cv.GPQR(Xc, pheno, group, Kfold = 5, taux = taux, check ="f1",method =penalty,plot.it=F)
  l=min(which(cvGS$cv==min(cvGS$cv)))
  nbeta <- c(cvGS$finalfit$b0[l], cvGS$finalfit$beta[,l])
  matplot(positionC/1e6,nbeta[-1],col=4, type = "p",pch=19,cex=1,ylab=ylab,xlab="Location",
          main= paste("GPQR -",penalty," tau=",taux),ylim=c(-0.008,0.006))
  abline(v=11.235,lty=2)
  abline(v=11.385,lty=2)
}
