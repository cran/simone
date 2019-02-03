rm(list=ls())
library(simone)
data(cancer)

res.coop <- simone(cancer$expr, tasks=cancer$status)
plot(res.coop, ask=FALSE)

glist <- getNetwork(res.coop, "BIC")
plot(glist[[1]],glist[[2]])
glist <- getNetwork(res.coop, 65)
plot(glist[[1]],glist[[2]])

