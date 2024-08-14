\name{VarScore}
\alias{VarScore}
\title{A score test for detecting extreme values in a vector autoregressive model}
\usage{
vech(x)
vec(x)
myVAR(samples, p=NULL)
Montecarlo(mu, sigma, B, sample_n)
gene_indep_variable(samples, p, k)
calc_I22_mean(z, k, p, n, sigma_inv)
calc_I22_var(z, k, p, n, sigma_inv)
mean_shift(B, k, p, samples, sigma, alpha)
variance_weight(B, k, p, samples, sigma, alpha)
}
\description{
A score test for detecting extreme values in a vector autoregressive model
}
\examples{
library(VarScore)
V1<-myVAR(IBM_SP)
new_sigma<-V1$hat_sigma
new_B<-V1$hat_B
residual<-data.frame(t(V1$hat_U))
k<-dim(new_B)[1]
p<-(dim(new_B)[2]-1)%/%k

S_mean_shift<-mean_shift(new_B, k, p, IBM_SP, new_sigma, 0.05)
ms_scores<-S_mean_shift$all
ms_outliters<-S_mean_shift$outliters
ms_correct<-S_mean_shift$correct

S_variance_weight<-variance_weight(new_B, k, p, IBM_SP, new_sigma, 0.05)
vw_scores<-S_variance_weight$all
vw_outliters<-S_variance_weight$outliters
vw_correct<-S_variance_weight$correct

par(mfrow=c(2,1))
n<-dim(samples)[1]
critical_ms<-qchisq(1-0.05/(n-p), k)
critical_vw<-qchisq(1-0.05/(n-p), 1)
plot(ms_scores$ms_score~as.Date(samples[1:dim(samples)[1],1]),type='l',col='blue',lwd=2,main='Score statistical value of MSOM',xlab='',ylab='')
abline(h=critical_ms,col='red',lwd=2)
legend("topright", legend=c("score","benchmark"), col=c('blue',"red"), lty=c(1,1), lwd=2)
plot(vw_scores$vw_score~as.Date(samples[1:dim(samples)[1],1]),type='l',col='blue',lwd=2,main='Score statistical value of CWM',xlab='',ylab='')
abline(h=critical_vw,col='red',lwd=2)
legend("topright", legend=c("score","benchmark"), col=c('blue',"red"), lty=c(1,1), lwd=2)

index<-c(row.names(ms_outliters),row.names(vw_outliters))
total_outliters_index<-as.numeric(index[!duplicated(index)])
result<-samples[sort(total_outliters_index),]
result$MSOM<-''
result$CWM<-''
result[row.names(ms_outliters),]$MSOM<-'*'
result[row.names(vw_outliters),]$CWM<-'*'
write.csv(result,"result.csv", row.names = FALSE)
}
