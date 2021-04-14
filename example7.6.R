# plot the log(Phi(x))

set.seed(197789);
# pdf("quardratic_approximation.pdf")


# par(mfrow=c(2,2));
# x=rnorm(10000);
# y1=log(pnorm(0.1*x));
# u=predict(lm(y1~x+I(x^2)));
# data=cbind(x,y1,u);
# data_sort=data[order(x),];

# matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
# xlab="X", ylab="log(Phi(k*X))");
# title(main="Quadratic approximation, k=0.1", cex=0.5);

# y1=log(pnorm(0.5*x));
# u=predict(lm(y1~x+I(x^2)));
# data=cbind(x,y1,u);
# data_sort=data[order(x),];
# matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
# xlab="X", ylab="log(Phi(k*X))", main="Quadratic approximation, k=0.5");

# y1=log(pnorm(1*x));
# u=predict(lm(y1~x+I(x^2)));
# data=cbind(x,y1,u);
# data_sort=data[order(x),];
# matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
# xlab="X", ylab="log(Phi(k*X))", main="Quadratic approximation, k=1");

# y1=log(pnorm(2*x));
# u=predict(lm(y1~x+I(x^2)));
# data=cbind(x,y1,u);
# data_sort=data[order(x),];
# matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
# xlab="X", ylab="log(Phi(k*X))", main="Quadratic approximation, k=2");


# dev.off();


# plot(x, y1, type="l");
# plot(x, u);

# matplot(x, cbind(y1,u));

# pdf("linear_approximation.pdf")

# par(mfrow=c(2,2));
# x=rnorm(10000);
# y2=log(pnorm(0.1*x+0.1)/pnorm(0.1*x));
# u=predict(lm(y2~x));
# data=cbind(x,y2,u);
# data_sort=data[order(x),];

# matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
# xlab="X", ylab="log(Phi(k*X+k)/Phi(k*X))");
# title(main="Linear approximation, k=0.1");

# y2=log(pnorm(0.5*x+0.5)/pnorm(0.5*x));
# u=predict(lm(y2~x));
# data=cbind(x,y2,u);
# data_sort=data[order(x),];

# matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
# xlab="X", ylab="log(Phi(k*X+k)/Phi(k*X))");
# title(main="Linear approximation, k=0.5");

# y2=log(pnorm(1*x+1)/pnorm(1*x));
# u=predict(lm(y2~x));
# data=cbind(x,y2,u);
# data_sort=data[order(x),];

# matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
# xlab="X", ylab="log(Phi(k*X+k)/Phi(k*X))");
# title(main="Linear approximation, k=1");

# y2=log(pnorm(2*x+2)/pnorm(2*x));
# u=predict(lm(y2~x));
# data=cbind(x,y2,u);
# data_sort=data[order(x),];

# matplot(x=data_sort[,1], y=data_sort[,2:3], type=rep("l",2), lwd=2, 
# xlab="X", ylab="log(Phi(k*X+k)/Phi(k*X))");
# title(main="Linear approximation, k=2");

# dev.off();

# pdf("quadratic_approximation_2.pdf")

par(mfrow=c(1,1));

x=rnorm(10000);
beta0=0.1;
beta1=0.1;
# mu=mean(x);
# beta0=0.5;
# beta1=0.5;

# beta0=1;
# beta1=1;

# beta0=2;
# beta1=2;

# beta0=-3;
# beta1=-3;

# plot for example 7.6
y1=log(1+exp(beta0+beta1*x));
# x_mean=mean(exp(beta0+beta1*x)/(1+exp(beta0+beta1*x)));

# u=log(1+exp(beta0+beta1*mu))+(x-mu)*beta1*exp(beta0+beta1*mu)/(1+exp(beta0+beta1*mu))+(x-mu)^2/2*beta1^2*exp(beta0+beta1*mu)/(1+exp(beta0+beta1*mu))^2;


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


