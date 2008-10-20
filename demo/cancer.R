library(simone)

## load the data set
data(cancer)

p <- length(gene.names)

cat("\n Breast cancer data set: real network with", p, "genes")

## SIMoNe parameters
Q <- 2
mult <- list(init=1,intra=1,inter=2,dust=4)
penalty.pcr <- Penalty(cancer.pcr, risk=0.1*p)[1,2]
penalty.not <- Penalty(cancer.notpcr, risk=0.1*p)[1,2]
rho.fracs <- c(1e-2,5e-2,7.5e-2,.25)
rhos.pcr  <- rho.fracs * penalty.pcr
rhos.not  <- rho.fracs * penalty.not

for (i in 1:length(rho.fracs)) {
  
  cat("\n Inferring network on the 'pCR' sample with rho =",rhos.pcr[i])
  res.pcr <- simone(cancer.pcr,Q,rho=rhos.pcr[i],multipliers=mult )
  cat("\n Inferring network on the 'not pCR' sample with rho =",rhos.not[i])
  res.not <- simone(cancer.notpcr,Q,rho=rhos.not[i],multipliers=mult)
  
  a <- readline(prompt="\nPress enter to continue...")
  
  layout(matrix(c(1,2),1,2),width=c(1,1))
  Gplot(res.pcr$K.hat, res.pcr$cl, cols=c(grey(0),grey(.8),grey(1)),
        display.labels=TRUE, labels=gene.names,main="pCR")
  Gplot(res.not$K.hat, res.not$cl, cols=c(grey(0),grey(.8),grey(1)),
        display.labels=TRUE,labels=gene.names,main="not pCR")
}
