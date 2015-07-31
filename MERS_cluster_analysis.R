############################################################################
#                                                                          #
#   THE ROLE OF SUPERSPREADING IN MIDDLE EAST RESPIRATORY                  #
#   SYNDROME CORONAVIRUS (MERS-COV) TRANSMISSION                           #
#                                                                          #
#   REFERENCE: KUCHARSKI AJ, ALTHAUS CL. EURO SURVEILL. 2015;20:pii=21167. #
#   WWW: http://www.eurosurveillance.org/ViewArticle.aspx?ArticleId=21167  #
#                                                                          #
#   Copyright (C) 2015 by Adam J. Kucharski (adam.kucharski@lshtm.ac.uk)   #
#   and Christian L. Althaus (christian.althaus@alumni.ethz.ch)            #
#                                                                          #
#   This program is free software; you can redistribute it and/or modify   #
#   it under the terms of the GNU General Public License as published by   #
#   the Free Software Foundation; either version 2 of the License, or      #
#   (at your option) any later version.                                    #
#                                                                          #
#   This program is distributed in the hope that it will be useful,        #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#   GNU General Public License for more details.                           #
#                                                                          #
#   You should have received a copy of the GNU General Public License      #
#   along with this program; if not, write to the                          #
#   Free Software Foundation, Inc.,                                        #
#   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              #
#                                                                          #
############################################################################

# Compile cluster data

# Cauchemez et al (2014)
chaindataMERS<-c(rep(1,27),rep(2,2),rep(3,4),rep(4,3),rep(5,2),7,13,26) # Table 1 scenario 2

# Breban et al (2013)
chaindataMERS_early<-c(rep(1,11),rep(2,2),rep(3,3),4,rep(5,2),24) # Table 1 scenario 2

# Poletto et al (2014)
chaindataMERS_2<-c(rep(1,42),rep(2,7),rep(3,2),5,5,10,22) # With Jordan cluster

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Functions to estimate R and k

# Outbreak size distribution function
hg.rj <- function(L,R,k) {
  p <- exp(lgamma(k*L+L-1) - lgamma(k*L) - lgamma(L+1) + log((R/k)^(L-1)/(1+R/k)^(k*L+L-1)))
  return(p)
}

# Likelihood function
likelihoodHG<- function(simdata,rstar,kstar){
  runs=length(simdata)
  liks=rep(1,runs)
  for(z in 1:runs){liks[z]=hg.rj(simdata[z],rstar,kstar)}
  sum(log(liks))
}

# Set up grid search
R0r1=0.25; R0r2=1.05; R0range=seq(R0r1,R0r2,by=.01)
kkr1=0.04; kkr2=55; kkrange=seq(kkr1,kkr2,by=0.01)

# Calculate likelihood surface
R0estimateHG<- function(simdata,R0range,kkrange){
  collect1 <- data.frame(matrix(NA, nrow=length(R0range),length(kkrange)))
  for(ii in 1:length(R0range)){
    for(jj in 1:length(kkrange)){
      rr=R0range[ii]
      kkr=kkrange[jj]
      collect1[ii,jj]=likelihoodHG(simdata,rr,kkr)
    }
  }
  collect1
}

liksurfMERS=R0estimateHG(chaindataMERS,R0range,kkrange)
liksurfMERS_2=R0estimateHG(chaindataMERS_2,R0range,kkrange)
liksurfMERSearly=R0estimateHG(chaindataMERS_early,R0range,kkrange)

likmMERS=(liksurfMERS==max(liksurfMERS))
likmMERS_2=(liksurfMERS_2==max(liksurfMERS_2))
likmMERSearly=(liksurfMERSearly==max(liksurfMERSearly))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate profile likelihoods

# For reference SARS 0.16 (90% confidence interval 0.11–0.64)
calc_profile<-function(liksurfMX,likmMX,chiV){
  
  prfk=apply(liksurfMX,2,function(x){max(x)})
  prfk2=kkrange[prfk-max(prfk)>-chiV]
  prfR=apply(liksurfMX,1,function(x){max(x)})
  prfR2=R0range[prfR-max(prfR)>-chiV]
  
  c(
    paste(R0range[sum(seq(1,length(R0range))%*%likmMX)]," (",min(prfR2),"-",max(prfR2),")",sep=""),
    paste(kkrange[sum(likmMX%*%seq(1,length(kkrange)))]," (",min(prfk2),"-",ifelse(max(prfk2)==max(kkrange),"Inf",max(prfk2)),")",sep="")
  )
  
}

output_parameters<-function(conf.interval){
  
  chiV=qchisq(conf.interval/100, df=1)/2 
  outputM1=calc_profile(liksurfMERS,likmMERS,chiV)
  outputME=calc_profile(liksurfMERSearly,likmMERSearly,chiV)
  outputM2=calc_profile(liksurfMERS_2,likmMERS_2,chiV)
  
  outputdata=data.frame(rbind(outputME,outputM1,outputM2))
  rownames(outputdata)=c("Breban","Cauchemez","Poletto")
  colnames(outputdata)=c("R","k")
  
  write.csv(outputdata,paste("param_estimates",conf.interval,".csv",sep=""))

}

output_parameters(90)
output_parameters(95)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot likelihood surface

colNN2=rgb(0.8,0.6,0)
colNN=rgb(0.1,0.4,0.8)
colNNg=rgb(0,0.7,0.1)

chiV90=qchisq(0.90, df=1)/2 
chiV95=qchisq(0.95, df=1)/2 

ctlnsMERS90<-contourLines(x=R0range,y=log10(kkrange),as.matrix(liksurfMERS-max(liksurfMERS)),levels=c(-chiV90,chiV90+0.001))
ctlnsMERS_290<-contourLines(x=R0range,y=log10(kkrange),as.matrix(liksurfMERS_2-max(liksurfMERS_2)),levels=c(-chiV90,chiV90+0.001))
ctlnsMERSearly90<-contourLines(x=R0range,y=log10(kkrange),as.matrix(liksurfMERSearly-max(liksurfMERSearly)),levels=c(-chiV90,chiV90+0.001))

ctlnsMERS<-contourLines(x=R0range,y=log10(kkrange),as.matrix(liksurfMERS-max(liksurfMERS)),levels=c(-chiV95,chiV95+0.001))
ctlnsMERS_2<-contourLines(x=R0range,y=log10(kkrange),as.matrix(liksurfMERS_2-max(liksurfMERS_2)),levels=c(-chiV95,chiV95+0.001))
ctlnsMERSearly<-contourLines(x=R0range,y=log10(kkrange),as.matrix(liksurfMERSearly-max(liksurfMERSearly)),levels=c(-chiV95,chiV95+0.001))

library(RColorBrewer)

ctlnsMERSstore<-contourLines(x=R0range,y=kkrange,as.matrix(liksurfMERS_2-max(liksurfMERS_2)),levels=c(-chiV95,chiV95+0.001))
write.csv(cbind(ctlnsMERSstore[[1]][[2]],ctlnsMERSstore[[1]][[3]]),"95_contour.csv")

xtick=seq(0,1.4,0.1)
ytick=c(0.01,0.05,0.2,0.1,0.5,1,2)
ytick=c(0.01,0.05,0.1,0.5,1,5,10,50,100)

filled.contour(x=R0range,y=log10(kkrange),as.matrix(liksurfMERS-max(liksurfMERS)),col='white',levels=seq(-3,0,0.5), #col=brewer.pal(8, "Blues")
               xlab=expression('Basic reproduction number, R'[0]),ylab="Dispersion parameter, k", 
               xlim=c(0.25,1.15),
               ylim=log10(c(0.05,max(kkrange))),
               plot.axes = {
                 axis(1, at=xtick, label=xtick);
                 axis(2, at=log10(ytick), label=ytick);
                 
                 lines(ctlnsMERSearly[[1]][[2]], ctlnsMERSearly[[1]][[3]], lwd=2,col=colNNg) 
                 lines(ctlnsMERSearly90[[1]][[2]], ctlnsMERSearly90[[1]][[3]], lwd=1,col=colNNg,lty=2) 

                 lines(ctlnsMERS[[1]][[2]], ctlnsMERS[[1]][[3]], lwd=2,col=colNN) 
                 lines(ctlnsMERS90[[1]][[2]], ctlnsMERS90[[1]][[3]], lwd=1,col=colNN,lty=2) 

                 lines(ctlnsMERS_2[[1]][[2]], ctlnsMERS_2[[1]][[3]], lwd=2,col=colNN2) 
                 lines(ctlnsMERS_290[[1]][[2]], ctlnsMERS_290[[1]][[3]], lwd=1,col=colNN2,lty=2) 
                 
                 points(R0range[sum(seq(1,length(R0range))%*%likmMERSearly)],log10(kkrange[sum(likmMERSearly%*%seq(1,length(kkrange)))]),pch=19,cex=1.5,col=colNNg)
                 points(R0range[sum(seq(1,length(R0range))%*%likmMERS)],log10(kkrange[sum(likmMERS%*%seq(1,length(kkrange)))]),pch=19,cex=1.5,col=colNN)
                 points(R0range[sum(seq(1,length(R0range))%*%likmMERS_2)],log10(kkrange[sum(likmMERS_2%*%seq(1,length(kkrange)))]),pch=19,cex=1.5,col=colNN2)
                 
                 text(x=1.03,y=log10(10),labels="To 21 Jun 2013",col=colNNg)
                 text(x=0.66,y=log10(1.5),labels="To 8 Aug 2013",col=colNN)
                 text(x=0.7,y=log10(0.075),labels="To 31 Aug 2013",col=colNN2)
                 lines(c(0,1.15),log10(c(0.16,0.16)),col='red',lwd=2,lty=2);
                 text(x=1,y=log10(0.2),labels="SARS",col='red')
               }
)

dev.copy(pdf,paste("Figure1.pdf",sep=""),width=10,height=8)
dev.off()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Contour lines of probability that a transmission cluster reaches at least 150 cases

# Read 95% confidence intervals for estimates of R0 and k (based on Poletto et al data)
cloud <- read.csv("95_contour.csv")

# Probability to observe a transmission chain of size L
prob <- function(L,R,k) {
  p <- exp(lgamma(k*L+L-1) - lgamma(k*L) - lgamma(L+1) + (L-1)*log(R/k)-(k*L+L-1)*log(1+R/k))
  return(p)
}

# Relationship between R0, k, and the probability to observe a cluster >= L
L <- 150
R0 <- seq(0.45,1.25,0.005)
k <- 10^seq(-1.3,0.3,0.005)
P <- matrix(NA,nrow=length(R0),ncol=length(k))

for(i in 1:length(R0)) {
	for(j in 1:length(k)) {
		P[i,j] <- 1 - sum(prob(1:(L-1),R0[i],k[j]))
	}
}

par(mfrow=c(1,2))
l <- c(-1e-15,1e-5,1e-3,5e-3,1e-2,2e-2,5e-2,1e-1,1.5e-1,1)
image(R0,k,P,xlab=bquote("Basic reproduction number" ~ italic("R")["0"]),ylab=bquote("Dispersion parameter" ~ italic("k")),log="y",col=rev(heat.colors(length(l)-1)),breaks=l,main="A                                                         ")
contour(R0,k,P,levels=l,labcex=0.8,add=TRUE)
points(0.47,0.26,pch=19)
lines(cloud[,2],cloud[,3],lwd=1,lty=2)

base <- P[,which(k==1)]

l <- c(-1e-15,0.5,1,1.2,3,1e1,1e2,1e4,1e20)
image(R0,k,P/base,xlab=bquote("Basic reproduction number" ~ italic("R")["0"]),ylab=bquote("Dispersion parameter" ~ italic("k")),log="y",col=rev(heat.colors(length(l)-1)),breaks=l,main="B                                                         ")
contour(R0,k,P/base,levels=l,labcex=0.8,add=TRUE)
points(0.47,0.26,pch=19)
lines(cloud[,2],cloud[,3],lwd=1,lty=2)

dev.copy(pdf,paste("Figure2.pdf",sep=""),width=9,height=5)
dev.off()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Integrate over likelihood surface to calculate probability of large outbreak

# Focus on MERS_2 - Use Poletto et al data
R0r1A=0.01; R0r2A=1.2; R0rangeA=seq(R0r1A,R0r2A,by=.01)
kkr1A=0.01; kkr2A=2; kkrangeA=seq(kkr1A,kkr2A,by=0.01)

# Calculate and normalise likelihood surface
liksurfMERS_2A=R0estimateHG(chaindataMERS_2,R0rangeA,kkrangeA)
normLik=liksurfMERS_2A-log(sum(exp(liksurfMERS_2A)))

# Integrate over parameter space to find expected probability greater than "more_than" given "intN" introductions
LargeProb<- function(normLik,R0rangeA,kkrangeA,intN,more_than){
  collectA <- data.frame(matrix(NA, nrow=length(R0rangeA),length(kkrangeA)))
  for(ii in 1:length(R0rangeA)){
    for(jj in 1:length(kkrangeA)){
      rr=R0rangeA[ii]
      kk=kkrangeA[jj]
      Pk=1-sum(prob(1:(more_than-1),rr,kk))
      collectA[ii,jj]=exp(normLik[ii,jj])*(1-(1-Pk)^intN)

    }
  }
  sum(collectA)
}

# Estimate probability of outbreak >=150 for different numbers of introductions
intro <- c(1e2,5e2,1e3,2e3)

sapply(intro,function(x){LargeProb(normLik,R0rangeA,kkrangeA,x,150)})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probability to observe a cluster >= L, using the likelihood surface for R0 and k, and N introductions

# Calculate and normalise surface for k=1
liksurfMERS_2A_k1=R0estimateHG(chaindataMERS_2,R0rangeA,1)
normLik_k1=liksurfMERS_2A_k1-log(sum(exp(liksurfMERS_2A_k1)))

par(mfrow=c(1,2))
maxlength <- 200
P <- 1
for(i in 2:maxlength) {
  P <- c(P,LargeProb(normLik_k1,R0rangeA,1,1,i))
}

Pk <- 1
for(i in 2:maxlength) {
  Pk <- c(Pk,LargeProb(normLik,R0rangeA,kkrangeA,intN=1,more_than=i))
}

plot(1:maxlength,P,ty="l",lty=2,frame=FALSE,log="y",xlab="Cluster size",ylab="Probability",ylim=c(1e-6,1),axes=FALSE,col="blue",main="A                                                         ")
axis(1,c(1,50,100,150,200),c(1,50,100,150,200))
axis(2,c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1),c(1e-6,NA,1e-4,NA,1e-2,NA,1))
lines(1:maxlength,Pk,col="red")
legend("topright",inset=0,legend=c("k = 1","k = 0.26 (95% CI: 0.09-1.24)"),col=c("blue","red"),lty=c(2,1),bty="n")

intro <- c(1e2,5e2,1e3,2e3)
plot(1:maxlength,rep(NA,maxlength),frame=FALSE,xlab="Cluster size",ylab="Probability",ylim=c(0,1),axes=FALSE,col="blue",main="B                                                         ")
axis(1,c(1,50,100,150,200),c(1,50,100,150,200))
axis(2,seq(0,1,0.2),seq(0,1,0.2))
for(jj in intro) {
  p_intro <- 1
  for(i in 2:maxlength) {
    p_intro <- c(p_intro,LargeProb(normLik,R0rangeA,kkrangeA,intN=jj,more_than=i))
  }
  lines(1:maxlength,p_intro,col="red",lty=which(intro==jj))
}
legend("topright",inset=0,legend=c("N = 100","N = 500","N = 1000","N = 2000"),col=rep("red",4),lty=1:4,bty="n")

dev.copy(pdf,paste("Figure3.pdf",sep=""),width=12,height=7)
dev.off()
