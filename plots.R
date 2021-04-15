library(coda)
library(lattice)
library(gridExtra)
library(MASS)
library(Matrix)
library(mice)
library(mitools)
library(survival)
# source('datagen1.R')


nm1nc1b <- function( n, beta, rho, a,d, c,tao )
{
  ok <- 0
  while( ok == 0 )
    {
### generating the data
### the covariates
zz=mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),2,2))
zm=zz[,1]
zc=zz[,2]
z <- cbind(zm, zc)
### 
risk <- exp( z %*% beta )
u1 <- runif(n)
t <- c * ( -log(1-u1) / risk )^(1/tao)
u2 <- runif(n)
ctime <- c * ( -log(1-u2) )

case <- rep(0, n)
case[ t <= ctime ] <- 1

x <- t
x[case==0] <- ctime[case==0]

temp <- quantile( x, 0.95 )
case[ x >= temp ] <- 0

###  assign the conditional selection probability
sele.prob <- 1/(1 + exp(a + d*zc ))
u1 <- runif( n )
sele.ind <- rep(0, n)
sele.ind[ u1 < sele.prob ] <- 1 
   ok <- ifelse( min(table(sele.ind, case)) >= 5, 1, 0 ) 
  }
data <- data.frame( x, zm, zc, case, sele.ind ) 

return( data )
}


# panel_fn <- function(x, y, ...)
# {
#    panel.xyplot(x, y, ...)
#    panel.xyplot(x, y, type="smooth", col="red", lwd=2, ...)
# }



set.seed(999)
n <- 1000
beta <- c( 1, 1 )
rho <- 0.5



data1 <- matrix( 0, nrow = n, ncol = 5 )

    data1 <- nm1nc1b(n, beta, rho, 0,-1, 1,2)
    t <- data1$x
    zm <- data1$zm
    zc <- data1$zc
    case <- data1$case
    sele.ind <- data1$sele.ind

Ht=nelsonaalen(data1, x, case)


plot(t, zm, xlab="t", ylab="X1");

test=lm(zm~t);
test2=lowess(zm ~t);
abline(test, col="green", lty=2, lwd=3);
# abline(test2, col="blue", lty=2, lwd=3);
lines(test2$x,test2$y,col="red",lwd=3)


scatter.smooth(t, zm, col="blue", lwd=3, xlab="t", ylab="X1");
abline(test, col="green",  lty=2, lwd=3);

log_t=log(t);
plot(log_t, zm, xlab="log(t)", ylab="X1");
test=lm(zm~log_t);
test2=lowess(zm ~log_t);
abline(test, col="green", lty=2, lwd=3);
lines(test2$x,test2$y,col="red",lwd=3)



# lines(test2$t,test2$zm,col="red",lwd=3)
scatter.smooth(log_t, zm, col="blue", lwd=3, xlab="log(t)", ylab="X1");
abline(test, col="green", lty=2, lwd=3);


plot(log(t), zm, xlab="log(t)");


plot(Ht, zm, xlab="H(t)", ylab="X1");

test=lm(zm~Ht);
test2=lowess(zm ~Ht);
abline(test, col="green", lty=2, lwd=3);
lines(test2$x,test2$y,col="red",lwd=3)


# lines(test2$t,test2$zm,col="red",lwd=3)
scatter.smooth(Ht, zm, col="blue", lwd=3, xlab="H(t)", ylab="X1");
abline(test, col="green", lty=2, lwd=3);



plot(Ht, zm, xlab="H(t)");

plot(t[case==1], zm[case==1]);
plot(log(t)[case==1], zm[case==1]);
plot(Ht[case==1], zm[case==1]);

plot(t[case==0], zm[case==0]);
plot(log(t)[case==0], zm[case==0]);
plot(Ht[case==0], zm[case==0]);



# plot1=xyplot(zm~t, panel=panel_fn, xlab="t", ylab="Zm", main="failure+censored", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(min(t),max(t)))
# plot2=xyplot(zm~log(t), panel=panel_fn, xlab="log(t)", ylab="Zm", main="failure+censored", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(min(log(t)),max(log(t))))
Ht=nelsonaalen(data1, x, case)
plot3=xyplot(zm~Ht, panel=panel_fn, xlab="H(t)", ylab="Zm", main="failure+censored", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(0,1))

plot4=xyplot(zm[case==1]~t[case==1], panel=panel_fn, xlab="t", ylab="Zm", main="failure", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(min(t),max(t)))
plot5=xyplot(zm[case==1]~log(t)[case==1], panel=panel_fn, xlab="log(t)", ylab="Zm", main="failure", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(min(log(t)),max(log(t))))
Ht=nelsonaalen(data1, x, case)
plot6=xyplot(zm[case==1]~Ht[case==1], panel=panel_fn, xlab="H(t)", ylab="Zm", main="failure", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(0,1))

plot7=xyplot(zm[case==0]~t[case==0], panel=panel_fn, xlab="t", ylab="Zm", main="censored", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(min(t),max(t)))
plot8=xyplot(zm[case==0]~log(t)[case==0], panel=panel_fn, xlab="log(t)", ylab="Zm", main="censored", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(min(log(t)),max(log(t))))
Ht=nelsonaalen(data1, x, case)
plot9=xyplot(zm[case==0]~Ht[case==0], panel=panel_fn, xlab="H(t)", ylab="Zm", main="censored", cex.lab=1.5, ylim=c(min(zm),max(zm)), xlim=c(0,1))
grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9, nrow=3, ncol=3)



