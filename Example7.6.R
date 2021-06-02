# Example 7.6
# plot the log(Phi(x))

set.seed(197789);

# Fig. 7.1

x=rnorm(10000);
beta0=0.1;
beta1=0.1;

y1=log(1+exp(beta0+beta1*x));

u=predict(lm(y1~x+I(x^2)));
data=cbind(x,y1,u);
data_sort=data[order(x),];

matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
xlab="X", ylab="log(1+exp(beta0+beta1*X))", main="beta0=0.1, beta1=0.1", cex=0.5);

beta0=-0.1;
beta1=0.1;
y1=log(1+exp(beta0+beta1*x));
u=predict(lm(y1~x+I(x^2)));
data=cbind(x,y1,u);
data_sort=data[order(x),];
matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
xlab="X", ylab="log(1+exp(beta0+beta1*X))", main="beta0=-0.1, beta1=0.1");

beta0=1;
beta1=1;
y1=log(1+exp(beta0+beta1*x));
u=predict(lm(y1~x+I(x^2)));
data=cbind(x,y1,u);
data_sort=data[order(x),];
matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
xlab="X", ylab="log(1+exp(beta0+beta1*X))", main="beta0=1, beta1=1");

beta0=-1;
beta1=1;
y1=log(1+exp(beta0+beta1*x));
u=predict(lm(y1~x+I(x^2)));
data=cbind(x,y1,u);
data_sort=data[order(x),];
matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
xlab="X", ylab="log(1+exp(beta0+beta1*X))", main="beta0=-1, beta1=1");

beta0=2;
beta1=2;
y1=log(1+exp(beta0+beta1*x));
u=predict(lm(y1~x+I(x^2)));
data=cbind(x,y1,u);
data_sort=data[order(x),];
matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
xlab="X", ylab="log(1+exp(beta0+beta1*X))", main="beta0=2, beta1=2");

beta0=-2;
beta1=2;
y1=log(1+exp(beta0+beta1*x));
u=predict(lm(y1~x+I(x^2)));
data=cbind(x,y1,u);
data_sort=data[order(x),];
matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
xlab="X", ylab="log(1+exp(beta0+beta1*X))", main="beta0=-2, beta1=2");

beta0=3;
beta1=3;
y1=log(1+exp(beta0+beta1*x));
u=predict(lm(y1~x+I(x^2)));
data=cbind(x,y1,u);
data_sort=data[order(x),];
matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
xlab="X", ylab="log(1+exp(beta0+beta1*X))", main="beta0=3, beta1=3");

beta0=-3;
beta1=3;
y1=log(1+exp(beta0+beta1*x));
u=predict(lm(y1~x+I(x^2)));
data=cbind(x,y1,u);
data_sort=data[order(x),];
matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
xlab="X", ylab="log(1+exp(beta0+beta1*X))", main="beta0=-3, beta1=3");


