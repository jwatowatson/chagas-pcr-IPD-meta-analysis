library(expint)

my_lambdas = seq(from=log10(1/5000), to=log10(5000/5000),length.out=100)
par(las=1)
plot(my_lambdas, 1-expint::gammainc(1,5*(10^my_lambdas))^15,
     type='l', panel.first=grid(), xaxt='n', xlab='Parasite density per ml',
     ylab='Probability at least one draw has at least 1 parasite',lwd=2)
vals = c(0.0002, 0.001, 0.01, .1, 1, 10)
axis(1, at = log10(vals),labels = vals)
abline(h=0.5,lty=2)

lines(my_lambdas, 1-expint::gammainc(1,5*(10^my_lambdas))^5,lty=2,lwd=2)

lines(my_lambdas, 1-expint::gammainc(1,5*(10^my_lambdas)),lty=3,lwd=2)
