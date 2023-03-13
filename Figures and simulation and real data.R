x <- seq(-6,6,length=400)
library(VGAM)
library(stats4)
library(splines)
library(gmp)


pdfmed <-  function(x,n,r,mu,b){
  con <- factorialZ(n)/(factorialZ(r-1)*factorialZ(n-r))
  
  val <- (con)*dlaplace(x,mu,b)*(plaplace(x,mu,b)^(r-1))*((1-plaplace(x,mu,b))^(n-r))
  value <- as.double(val)
}



y <- pdfmed(x,171,(171-1)/2+1,4,7)

plot(x,y,type='l')


cdfmed <- function(n,r,mu,b){
  val <- integrate(Vectorize(pdfmed),lower = -100,upper = 100,n=n,r=r,mu=mu,b=b)$value
}

v <- cdfmed(171,86,4,7)
n <- 401
r <- (n-1)/2+1

pdf_l <-  function(x){
  con <- factorialZ(n)/(factorialZ(r-1)*factorialZ(n-r))
  
  val <- (con)*dlaplace(x,0,1)*(plaplace(x,0,1)^(r-1))*((1-plaplace(x,0,1))^(n-r))
  value <- as.double(val)
}
pdfz <- function(z,n,mu,b){
  r <- (n-1)/2+1
  val <- pdfmed(z*b+mu,n,r,mu,b)*b
  return(val)
}

z <- seq(-1,1,length=400)
w <- pdfz(z,n,4,7)
j <- pdf_l(z)
plot(z,w,type='l')
lines(z,j,type='l',col=2)

cdfz <- function(n,mu,b){
  val <- integrate(Vectorize(pdfz),lower = -Inf,upper = Inf,n=n,mu=mu,b=b)$value
}

k <- cdfz(n,0,1)

################################# Figure 1 #######################


pdf_l <-  function(x,n){
  m<-(n+1)/2
  con <- factorialZ(n)/(factorialZ(m-1)*factorialZ(n-m))
  
  val <- (con)*dlaplace(x,0,1)*(plaplace(x,0,1)^(m-1))*((1-plaplace(x,0,1))^(m-1))
  value <- as.double(val)
}
w <- seq(-1,1,length=400)
rez <- pdf_l(z,15)

plot(w,rez,type='l',ylim=c(0,6),col=1,lwd=1.5,lty=1,ylab=expression(f[W]))
lines(w,pdf_l(z,35),type='l',col=2,lwd=1.5,lty=2)
lines(w,pdf_l(z,105),type='l',col=5,lwd=1.5,lty=3)
lines(w,pdf_l(z,205),type='l',col=6,lwd=1.5,lty=4)
legend(0.6, 5.5, legend=c("n=15", "n=35", "n=105", "n=205"),
       col=c("1", "2", "5", "6"), lty=1:4, cex=0.8,
       box.lty=0)








cdfl <- function(n){
  val <- integrate(Vectorize(pdf_l),lower = -Inf,upper = Inf,n=n)$value
}
kk <- cdfl(15)




################## Figure 2 ####################################
library(gmp)

w <- seq(0,150,length=1000)

pdf_W <- function(w,n){
  val <- dgamma(w^2,shape = n,rate = 1/2)*2*w
}



plot(w,pdf_W(w,300),type = 'l',col=1)


check.pdf_W <- function(n){
  val <- integrate(Vectorize(pdf_W),lower = 0,upper = 100,n=n)
}
cc <- check.pdf_W(300)
cc
pdf_U <-  function(u,n){
  valfunc <- function(w,u,n){
 
   inte <- dnorm(u/w,0,1)*pdf_W(w,n)/w
  }
  value <- integrate(Vectorize(valfunc),lower = 0,upper = 50,u=u,n=n)$value
  final_value <- value
  return(final_value)
}

pdf_baru <- function(u,n){
  val <- n*pdf_U(n*u,n)
}

tt<- pdf_baru(2,30)


y <- seq(-1,1,length=1000)
yy <- rep(0,1000)
for (i in 1:1000){
  yy[i]<-pdf_baru(y[i],30)
}
y1 <- rep(0,1000)
for (i in 1:1000){
  y1[i]<-pdf_baru(y[i],80)
}
y2 <- rep(0,1000)
for (i in 1:1000){
  y2[i]<-pdf_baru(y[i],150)
}
y3 <- rep(0,1000)
for (i in 1:1000){
  y3[i]<-pdf_baru(y[i],300)
}
plot(y,yy,type = 'l',col=1,lwd=1.5,lty=1,ylab=expression(f[bar(Y)]),ylim=c(0,5.5))
lines(y,y1,type='l',col=4,lwd=1.5,lty=2)
lines(y,y2,type='l',col=5,lwd=1.5,lty=3)
lines(y,y3,type='l',col=8,lwd=1.5,lty=4)
legend(0.6,5, legend=c("n=30", "n=80", "n=150", "n=300"),
       col=c("1", "4", "5", "8"), lty=1:4, cex=0.8,
       box.lty=0)





check.pdf_U <- function(n){
  val <- 2*integrate(Vectorize(pdf_baru),lower = 0,upper = 1000,n=n)$value
}
aa <- check.pdf_U(n)
aa

############################## Figure 4 ############################


n <- 20
pdfcheck <- function(x,n){
  val <- pdf_U(sqrt(2*n)*x,n)*sqrt(2*n)
}
y <- seq(-5,5,length=400)

z <- rep(0,400)
for (i in 1:400){
  z[i]<-pdfcheck(y[i],n)
}

n <- 100
z1 <- rep(0,400)
for (i in 1:400){
  z1[i]<-pdfcheck(y[i],n)
}

plot(y,dnorm(y,0,1),type = 'l',col=1,xlim=c(-5,5),ylim=c(0,0.45),ylab="density curves",lty=1)
lines(y,z,type = 'l',col=2,lty=2)
lines(y,z1,type = 'l',col=4,lty=3)
legend(1.3, 0.4, legend=c("Z", "sqrt(n/2)Bar_Y with n=20","sqrt(n/2)Bar_Y with n=100"),
       col=c("1", "2","4"), lty=1:3, cex=0.8,
       box.lty=0)



pdfhh <- function(u,n){
  val <- dnorm(u/sqrt(2*n),0,1)/sqrt(2*n)
}
n <- 10
y <- rep(0,10000)
for (i in 1:10000){
  y[i]<-pdf_U(u[i],n)
}
w <- rep(0,10000)
for (i in 1:10000){
  w[i]<-pdfhh(u[i],n)
}
plot(u,y,type = 'l',col=1,xlim=c(-50,50))
lines(u,w,type = 'l',col=2)



################################# coverage rate ########################

library("dplyr")

cr<- function(n, mu,b,f){
  M <- 10000
  estitux <- rep(0,M)
  for (i in 1:M){
    sam <- rlaplace(n,mu,b)
    estitux[i] <- median(sam)
  }
  es <- mean(estitux)
  sd <- sd(estitux)
  interlower <- estitux-sqrt(2)*b*f
  interupper <- estitux+sqrt(2)*b*f
  cr <- sum(mu>= interlower& mu<=interupper)/M

  
  return(c(cr,es,sd))
}
cr(219,5,2,0.1)
cr(155,5,2,0.1)
cr(103,5,2,0.15)
cr(73,5,2,0.15)
cr(61,5,2,0.2)
cr(43,5,2,0.2)
cr(41,5,2,0.25)
cr(29,5,2,0.25)
cr(31,5,2,0.3)
cr(21,5,2,0.3)




################################ real date example #########################
library(quantmod)
GSPC <- getSymbols("^GSPC", from = "2019-12-31", to = "2022-12-31", auto.assign = FALSE)
prices <- Cl(GSPC)
dim(prices)
plot(log(prices), col = 'blue', lwd = 2, ylab = "log-price", main = "S&P 500 index")
R_daily <- diff(log(prices))[-1]
plot(R_daily, col = 'blue', lwd = 1, ylab = "log-return", main = "S&P 500 index")

h <- hist(as.vector(R_daily), breaks = 100, prob = TRUE, col = "lightgray", 
          xlab = "return", main = "Histogram of log-returns for 3 years data",ylim = c(0,50))
xfit <- seq(min(R_daily), max(R_daily), length = 100) 
yfit <- dlaplace(xfit, location = median(R_daily), scale = sum(abs(coredata(R_daily)-median(R_daily)))/dim(R_daily)[1])
zfit <- dnorm(xfit,mean=mean(R_daily),sd <- sd(R_daily))
lines(xfit, yfit, col = "blue", lty=1,lwd=2)
lines(xfit, zfit, col = "green", lty=1,lwd=2)
legend(0.05, 40, legend=c("Laplace", "Normal"),
       col=c("blue", "green"), lty=1:1, cex=0.8,lwd=2,
       box.lty=0)



library(moments)
kurtosis(R_daily)
s <- sort(coredata(R_daily))
kurtosis(s[3:750])



sample <- sample(R_daily,219)

plot(sample, col = 'blue', lwd = 1, ylab = "log-return", main = "S&P 500 index")

h <- hist(as.vector(sample), breaks = 40, prob = TRUE, col = "lightgray", 
          xlab = "return", main = "Histogram of log-returns for the sample",ylim = c(0,50))
xfit <- seq(min(sample), max(sample), length = 100) 
yfit <- dlaplace(xfit, location = median(sample), scale = sum(abs(coredata(sample)-median(sample)))/dim(sample)[1])
zfit <- dnorm(xfit,mean=mean(sample),sd <- sd(sample))
lines(xfit, yfit, col = "blue", lty=1,lwd=2)
lines(xfit, zfit, col = "green", lty=1,lwd=2)
legend(0.05, 40, legend=c("Laplace", "Normal"),
       col=c("blue", "green"), lty=1:1, cex=0.8,lwd=2,
       box.lty=0)

