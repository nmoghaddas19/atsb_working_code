u = 0.09
b = 0.09
d = seq(0,1,0.001)

feed_rate <- (-u*(d-1)-sqrt(u^2*(d-1)^2-4*(d-1)*d*b*u))/(2*(d-1))
feed_rate <- f^2/(b*u) + f/b

plot(d, feed_rate, type="l", frame.plot = F, xlim=c(0,0.5), ylim=c(0,1))

u <- seq(0,1,0.01)
d <- 0.1
b <- 0.3
f <- u*(d+u)/(b*d+u*(d+u))
plot(f, u, type="l", frame.plot = F, xlim=c(0,0.6), ylim=c(0,0.5), lwd=2,
     xlab="Dyed fraction", ylab="Feeding rate")
approx(f, u, xout = 0.35)

u <- seq(0,1,0.01)
r <- 1/4.5
d <- 0.1
b <- 0.1
f <- u*(d+u)/((b+r)*(d+r)+u*(d+u))
plot(f, u, type="l", frame.plot = F, xlim=c(0,0.6), ylim=c(0,0.5), lwd=2,
     xlab="Dyed fraction", ylab="Feeding rate")

a <- -1.595
c <- 1.824
f <- a*u^2 + c*u
lines(f, u, type="l", xlim=c(0,0.6), ylim=c(0,0.5), lwd=2,
     xlab="Dyed fraction", ylab="Feeding rate", col=2)

approx(f, u, xout = 0.35)
