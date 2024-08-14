library(MASS)
library(vars)
library(matrixcalc)

#' Return vech
#'
#' @param x the square matrix
#'
#' @return vech
#' @export
#'
#' @examples
#' library(VarScore)
#' a <- matrix(1:9,3,3)
#' vech(a)
vech<-function (x){
  t(x[lower.tri(x,diag=T)])
}

#' Return vec
#'
#' @param x the square matrix
#'
#' @return vec
#' @export
#'
#' @examples
#' library(VarScore)
#' a <- matrix(1:9,3,3)
#' vech(a)
vec<-function (x){
  t(t(as.vector(x)))
}

#' The ML estimators of B and Σ
#'
#' @param samples should be a data frame. Must be greater than three columns, and the first column is a timestamp.
#' @param p lags. default 'NULL'.
#'
#' @return The list include the ML estimators of B and Σ.
#' @export
#'
#' @examples
#' library(VarScore)
#' V1<-myVAR(IBM_SP)
#' new_sigma<-V1$hat_sigma
#' new_B<-V1$hat_B
#' residual<-data.frame(t(V1$hat_U))
myVAR<-function(samples, p=NULL){
  sample<-as.matrix(samples[2:dim(samples)[2]])
  if (length(p)==0){
    p<-VARselect(sample,type = "const")$selection["AIC(n)"]
  }
  k<-dim(sample)[2]
  n<-dim(sample)[1]
  Y<-t(sample[n:(p+1),])
  X1<-matrix(1,1,n-p)
  X2<-matrix(0,k*p,n-p)
  for (i in 1:(n-p)) {
    X2[,i]<-vec(t(matrix(array(sample[(n-i):(n-p+1-i),]),p,k)))
  }
  X<-rbind(X1,X2)
  hat_B<-Y%*%t(X)%*%solve(X%*%t(X))
  hat_U<-Y-hat_B%*%X
  hat_sigma<-1/(n-p)*hat_U%*%t(hat_U)
  return(list(hat_B=hat_B,hat_U=hat_U,hat_sigma=hat_sigma))
}

#' Montecarlo
#'
#' @param mu The mean of multivariate normal distribution.
#' @param sigma The Σ of multivariate normal distribution.
#' @param B The coefficient of VAR.
#' @param sample_n The number of simulation.
#'
#' @return data frame
#' @export
#'
#' @examples
#' library(VarScore)
#' mu <- matrix(c(0.003,0.006,-0.015),3,1)
#' sigma <- matrix(c(1,0.6,0.76,0.6,1,0.14,0.76,0.14,1),3,3)
#' B <- matrix(c(-0.338,-0.134,-2.221,0.746,0.328,4.586,0.058,0.043,0.301),3,3)
Montecarlo <- function(mu, sigma, B, sample_n){
  B <- cbind(mu,B)
  vec_dimension <- length(mu)
  l <- matrix(1,1,1)

  mrn <- t(mvrnorm(sample_n,matrix(0,vec_dimension,1),sigma))

  p <- (dim(B)[2]-1)%/%vec_dimension
  x <- matrix(rep(0,sample_n*vec_dimension),sample_n,vec_dimension)
  x_patial <- rbind(l,t(t(array(t(x[(sample_n-p+1):sample_n,])))))

  for (i in 1:(sample_n-p)) {
    new_x <- B%*%x_patial+mrn[,sample_n-p-i+1]
    x[(sample_n-p-i+1),]<-new_x
    x_patial[2:dim(x_patial)[1]]<-array(t(x[(sample_n-p-i+1):(sample_n-i),]))
  }
  colnames(x)<-paste("seq",1:vec_dimension)
  data<-data.frame(date=1:sample_n,x[,1:vec_dimension][sample_n:1,])
  return(data)
}

#' The matrix of independent variables
#'
#' @param samples should be a data frame. Must be greater than three columns, and the first column is a timestamp.
#' @param p Lags.
#' @param k The dimension of VAR.
#'
#' @return The matrix of independent variables
#' @export
#'
#' @examples
#' library(VarScore)
#' gene_indep_variable(IBM_SP, 1, 2)
gene_indep_variable <- function(samples, p, k){
  sample<-as.matrix(samples[2:dim(samples)[2]])
  n <- dim(sample)[1]
  l <- matrix(1,1,1)
  z <- matrix(0,n-p,k*p+1)
  for (i in 1:(n-p)) {
    temp <- cbind(l,t(array(t(sample[,1:dim(sample)[2]][(i+p-1):i,]))))
    z[i,]<-temp
  }
  return(z)
}

#' The I22 of Mean-shift perturbation model
#'
#' @param z The matrix of independent variables.
#' @param k The dimension of VAR.
#' @param p Lags.
#' @param n The number of samples.
#' @param sigma_inv The inverse matrix of Σ.
#'
#' @return The I22 of Mean-shift perturbation model
#' @export
#'
#' @examples
#' library(VarScore)
#' x<-gene_indep_variable(IBM_SP, 1, 2)
#' V<-myVAR(IBM_SP)
#' sigma<-V$hat_sigma
#' B<-V$hat_B
#' sigma_inv<-solve(sigma)
#' k<-dim(sigma)[1]
#' p<-(dim(B)[2]-1)%/%k
#' n<-dim(IBM_SP)[1]
#' calc_I22_mean(x, k, p, n, sigma_inv)
calc_I22_mean <- function(z, k, p, n, sigma_inv){
  dup_mat <- D.matrix(k)
  I22 <- matrix(0,((1+2*p)*k*k+3*k)%/%2,((1+2*p)*k*k+3*k)%/%2)
  I221 <- kronecker(t(z)%*%z, sigma_inv)
  I222 <- (n-p)/2*t(dup_mat)%*%kronecker(sigma_inv,sigma_inv)%*%dup_mat
  I22[1:(k*k*p+k),1:(k*k*p+k)]<-I221
  I22[(k*k*p+k+1):dim(I22)[1],(k*k*p+k+1):dim(I22)[2]]<-I222
  return(I22)
}

#' The I22 of Case-weight perturbation model
#'
#' @param z The matrix of independent variables.
#' @param k The dimension of VAR.
#' @param p Lags.
#' @param n The number of samples.
#' @param sigma_inv The inverse matrix of Σ.
#'
#' @return The I22 of Case-weight perturbation model
#' @export
#'
#' @examples
#' library(VarScore)
#' x<-gene_indep_variable(IBM_SP, 1, 2)
#' V<-myVAR(IBM_SP)
#' sigma<-V$hat_sigma
#' B<-V$hat_B
#' sigma_inv<-solve(sigma)
#' k<-dim(sigma)[1]
#' p<-(dim(B)[2]-1)%/%k
#' n<-dim(IBM_SP)[1]
#' calc_I22_var(x, k, p, n, sigma_inv)
calc_I22_var <- function(z, k, p, n, sigma_inv){
  dup_mat <- D.matrix(k)
  I22 <- matrix(0,((1+2*p)*k*k+3*k)%/%2,((1+2*p)*k*k+3*k)%/%2)
  I221 <- kronecker(t(z)%*%z, sigma_inv)
  I222 <- (n-p-1)/2*t(dup_mat)%*%kronecker(sigma_inv,sigma_inv)%*%dup_mat
  I22[1:(k*k*p+k),1:(k*k*p+k)]<-I221
  I22[(k*k*p+k+1):dim(I22)[1],(k*k*p+k+1):dim(I22)[2]]<-I222
  return(I22)
}

#' The Mean-shift perturbation model
#'
#' @param B The hat-B of VAR.
#' @param k The dimension of VAR.
#' @param p Lags.
#' @param samples should be a data frame. Must be greater than three columns, and the first column is a timestamp.
#' @param sigma The hat-sigma of VAR.
#' @param alpha The significance level.
#'
#' @return The result of Mean-shift perturbation model
#' @export
#'
#' @examples
#' library(VarScore)
#' V1<-myVAR(IBM_SP)
#' new_sigma<-V1$hat_sigma
#' new_B<-V1$hat_B
#' residual<-data.frame(t(V1$hat_U))
#' k<-dim(new_B)[1]
#' p<-(dim(new_B)[2]-1)%/%k
#'
#' S_mean_shift<-mean_shift(new_B, k, p, samples, new_sigma, 0.05)
#' ms_scores<-S_mean_shift$all
#' ms_outliters<-S_mean_shift$outliters
#' ms_correct<-S_mean_shift$correct
mean_shift<-function(B, k, p, samples, sigma, alpha){
  sigma_inv<-solve(sigma)
  z<-gene_indep_variable(samples, p, k)
  dup_mat <- D.matrix(k)
  sample<-as.matrix(samples[,2:dim(samples)[2]])
  n<-dim(sample)[1]
  I22<-calc_I22_mean(z, k, p, n, sigma_inv)
  critical<-qchisq(1-alpha/(n-p), k)
  score<-matrix(0,n-p,1)
  for (i in 1:(n-p)) {
    zi<-t(t(z[i,]))
    Lr <- sigma_inv%*%(t(t(sample[i+p,]))-B%*%zi)

    I12<-matrix(0,k,((1+2*p)*k*k+3*k)%/%2)
    I12[,1:(k*k*p+k)]<-sigma_inv%*%kronecker(t(zi),diag(k))

    I21<-t(I12)

    score[i]<-t(Lr)%*%solve(sigma_inv-I12%*%solve(I22)%*%I21)%*%Lr
  }
  samples$ms_score<-rbind(matrix(NA,p,1),score)
  return(list(all=samples,outliters=samples[which(samples$ms_score>=critical),],correct=samples[which(samples$ms_score<critical),]))
}

#' The Case-weight perturbation model
#'
#' @param B The hat-B of VAR.
#' @param k The dimension of VAR.
#' @param p Lags.
#' @param samples should be a data frame. Must be greater than three columns, and the first column is a timestamp.
#' @param sigma The hat-sigma of VAR.
#' @param alpha The significance level.
#'
#' @return The result of Case-weight perturbation model
#' @export
#'
#' @examples
#' library(VarScore)
#' V1<-myVAR(IBM_SP)
#' new_sigma<-V1$hat_sigma
#' new_B<-V1$hat_B
#' residual<-data.frame(t(V1$hat_U))
#' k<-dim(new_B)[1]
#' p<-(dim(new_B)[2]-1)%/%k
#'
#' S_variance_weight<-variance_weight(new_B, k, p, IBM_SP, new_sigma, 0.05)
#' vw_scores<-S_variance_weight$all
#' vw_outliters<-S_variance_weight$outliters
#' vw_correct<-S_variance_weight$correct
variance_weight<-function(B, k, p, samples, sigma, alpha){
  sigma_inv<-solve(sigma)
  z<-gene_indep_variable(samples, p, k)
  dup_mat <- D.matrix(k)
  sample<-as.matrix(samples[,2:dim(samples)[2]])
  n<-dim(sample)[1]
  I22<-calc_I22_var(z, k, p, n, sigma_inv)
  critical<-qchisq(1-alpha/(n-p), 1)
  score<-matrix(0,n-p,1)
  vec_sigma_inv<-vec(sigma_inv)

  I12 <- matrix(0,1,((1+2*p)*k*k+3*k)%/%2)
  I12[(k*k*p+k+1):dim(I12)[2]] <- 1/2*t(vec_sigma_inv)%*%dup_mat

  I21<-t(I12)

  b <- I12%*%solve(I22)%*%I21

  for (i in 1:(n-p)) {
    zi<-t(t(z[i,]))
    Lw_square<-((k-t(t(t(sample[i+p,]))-B%*%zi)%*%sigma_inv%*%(t(t(sample[i+p,]))-B%*%zi))/2)^2

    score[i]<-Lw_square*solve(k/2-b)
  }
  samples$vw_score<-rbind(matrix(NA,p,1),score)
  return(list(all=samples,outliters=samples[which(samples$vw_score>=critical),],correct=samples[which(samples$vw_score<critical),]))
}

