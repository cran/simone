library(simone)

## Data set and grapĥ generation
p <- 200
n <- 800
proba.in  <- 0.125
proba.out <- 0.00125
alpha <- c(1/3,1/3,1/3)
Q <- length(alpha)
X <- SimDataAffiliation (p, n, proba.in, proba.out, alpha)

## Try an adapted value for the base Penalty
expectedNumberOfEdges <- sum(abs(X$K.theo)>0)
rho.base <- Penalty(X$data,risk=0.25*expectedNumberOfEdges)[1,2]

## Compare SIMoNe to direct GLasso
res.Simone <- simone(X$data, Q, rho=rho.base, cl.theo=X$cl.theo, verbose=TRUE)
res.GLasso <- InferEdges(X$data,rho.base, method="glasso")

par(mfrow=c(2,2))
g <- Gplot(X$K.theo, X$cl.theo, main="Theoretical graph")
Gplot(res.Simone$K.hat, res.Simone$cl, coord=g, main="SIMoNe Inference")
Gplot(res.GLasso$K.hat, coord=g, main="GLasso Inference")
Gplot(res.Simone$K.hat.perfect, X$cl.theo, coord=g, main="Perfect SIMoNe")
