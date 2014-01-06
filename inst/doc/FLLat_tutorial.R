### R code from vignette source 'FLLat_tutorial.rnw'

###################################################
### code chunk number 1: FLLat_tutorial.rnw:105-109
###################################################
library(FLLat)
data(simaCGH)
result.pve <- FLLat.PVE(simaCGH,J.seq=1:ncol(simaCGH))
plot(result.pve)


###################################################
### code chunk number 2: FLLat_tutorial.rnw:131-134
###################################################
result.bic <- FLLat.BIC(simaCGH,J=5)
result.bic$lam1
result.bic$lam2


###################################################
### code chunk number 3: FLLat_tutorial.rnw:145-146
###################################################
plot(result.bic$opt.FLLat)


###################################################
### code chunk number 4: FLLat_tutorial.rnw:156-157
###################################################
plot(result.bic$opt.FLLat,what="weights")


###################################################
### code chunk number 5: FLLat_tutorial.rnw:217-220
###################################################
result.fdr <- FLLat.FDR(simaCGH,result.bic$opt.FLLat)
result.fdr$thresh.control
plot(result.fdr)


