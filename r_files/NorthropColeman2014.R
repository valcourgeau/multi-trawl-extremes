#########################################################################
# R code to accompany                                                   #
# Northrop, P.J. and Coleman, C.L.(2014) Improved threshold diagnostic  #
# plots for extreme value analyses, Extremes, To appear.                #
#########################################################################

# Set working directory 
# Should contain NorthropColeman2014.fns and NJ2011GoMdata.txt.
setwd("C:/Users/paulnorthrop/Documents/UGPROJ/CLAIRE/PAPER R/SWEAVE")

# Reads in functions and nidd.thresh data
source("NorthropColeman2014.fns")

#-------------------------- Figure 4 -------------------------------#

u <- c(seq(from=65,to=115,by=5),120)
# Score test only
res.nidd.5 <- score.fitrange(nidd.thresh,u)               
# score and LR tests (takes much longer)
res.nidd.5.LRT <- score.fitrange(nidd.thresh,u,do.LRT=T)  

# Look at output object 
# [$nllh is the negated log-likelihood at the (restricted) MLE.]
res.nidd.5.LRT

# Figure 4a (excluding lines based on simulation) ...
 par(mar=c(4.2,4.2,2.2,1),lwd=2,cex.lab=1.5,cex.axis=1.5)
 my.xlab <- expression(paste("lowest threshold / ", m^3*s^-1))
 plot(res.nidd.5.LRT$u,res.nidd.5.LRT$e.p.values,type="b",lty=1,pch="S",ylim=c(0,1),xlab=my.xlab,ylab="p-value")
 lines(res.nidd.5.LRT$u,res.nidd.5.LRT$LRT.p.values,type="b",lty=1,pch="L")
 lines(res.nidd.5.LRT$u,res.nidd.5.LRT$e.mult.mult,type="l",lty=2,pch="")
 lines(res.nidd.5.LRT$u,res.nidd.5.LRT$l.mult.mult,type="l",lty=4,pch="")
 axis(3,at=res.nidd.5.LRT$u,labels=res.nidd.5.LRT$n.between[1:length(res.nidd.5.LRT$u)],cex.axis=1)

# Score test only, fewer thresholds
u <- c(seq(from=65,to=115,by=10),120)
res.nidd.10 <- score.fitrange(nidd.thresh,u)

# Score test only, more thresholds
u <- c(seq(from=65,to=115,by=2.5),120)
res.nidd.2.5 <- score.fitrange(nidd.thresh,u)

# Figure 4b (excluding lines based on simulation) ...

pdf(file="figure4b.pdf",width=6*sqrt(2),height=6)
 par(mar=c(4.2,4.2,2.2,1),lwd=2,cex.lab=1.5,cex.axis=1.5)
 my.xlab <- expression(paste("lowest threshold / ", m^3*s^-1))
 plot(res.nidd.5$u,res.nidd.5$e.p.values,type="b",ylim=c(0,1),xlab=my.xlab,ylab="p-value",pch=16,lty=1)
 lines(res.nidd.10$u,res.nidd.10$e.p.values,type="b",pch=4,lty=2)
 lines(res.nidd.2.5$u,res.nidd.2.5$e.p.values,type="b",pch=1,lty=4)
 legend("topleft",legend=c("m=7","m=12","m=22"),lty=c(2,1,4),pch=c(4,16,1),cex=1.5)
dev.off()


