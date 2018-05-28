#p.adjust
x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
p <- 2*pnorm(sort(-abs(x)))

round(p, 3)
round(p.adjust(p), 3)
round(p.adjust(p, "BH"), 3)

## or all of them at once (dropping the "fdr" alias):
p.adjust.M <- p.adjust.methods[p.adjust.methods != "fdr"]
p.adj    <- sapply(p.adjust.M, function(meth) p.adjust(p, meth))
p.adj.60 <- sapply(p.adjust.M, function(meth) p.adjust(p, meth, n = 60))
stopifnot(identical(p.adj[,"none"], p), p.adj <= p.adj.60)
round(p.adj, 3)
## or a bit nicer:
noquote(apply(p.adj, 2, format.pval, digits = 3))


## and a graphic:
matplot(p, p.adj, ylab="p.adjust(p, meth)", type = "l", asp = 1, lty = 1:6,
        main = "P-value adjustments")
legend(0.7, 0.6, p.adjust.M, col = 1:6, lty = 1:6)

## Can work with NA's:
pN <- p; iN <- c(46, 47); pN[iN] <- NA
pN.a <- sapply(p.adjust.M, function(meth) p.adjust(pN, meth))
## The smallest 20 P-values all affected by the NA's :
round((pN.a / p.adj), 4)




#qvalue
# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")

library(qvalue)
data(hedenfalk)
p <- hedenfalk$p

# get q-value object
qobj <- qvalue(p)
plot(qobj)
hist(qobj)

# options available
qobj <- qvalue(p, lambda=0.5, pfdr=TRUE)
qobj <- qvalue(p, fdr.level=0.05, pi0.method="bootstrap", adj=1.2)


##multtest
# source("https://bioconductor.org/biocLite.R")
# biocLite("multtest")
library(multtest)
procedures = c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", 
               "BY")
adjusted = mt.rawp2adjp(p, procedures)
adjusted$adjp[1:10, ]
adjusted$adj[order(adjusted$index)[1:10], ]
p.adj    <- sapply(p.adjust.M, function(meth) p.adjust(p, meth))