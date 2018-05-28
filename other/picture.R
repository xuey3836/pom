
##boxplot
setwd("C:/Users/o0/Desktop/MSE for POM/code")
# pm=read.csv(file = "pm-no.csv",header = T)
pm=read.csv(file = "pm-with.csv",header = T)
par(mfrow=c(1,1))
beta=0
for (k in 1:3){
  boxplot(pm[,2:3+2*(k-1)], at = 1:2,col=c("3","4"),boxwex = 0.6,ylim=c(-0.4,0.4),xlim=c(0,9),xaxt = "n")
  abline(h=beta,lty = 2)
  boxplot(pm[,4:5+2*(k-1)], at = 1:2+3 ,col=c("3","4"),add=T,boxwex = 0.6,xaxt = "n")
  boxplot(pm[,6:7+2*(k-1)], at = 1:2+6,col=c("3","4"),add=T,boxwex = 0.6,xaxt = "n")
  axis(1, at = seq(1.5,9,3), labels = c("MAF=0.05","MAF=0.10","MAF=0.15"), tick = TRUE)
  axis(1, at = seq(1.5,9,3), labels = c("MAF=0.20","MAF=0.25","MAF=0.30"), tick = TRUE)
  axis(1, at = seq(1.5,9,3), labels = c("MAF=0.35","MAF=0.40","MAF=0.45"), tick = TRUE)
}

pm=read.csv(file = "pm-gamma.csv",header = T)
gamma=5
for (k in 1:3){
  boxplot(pm[,2:3+2*(k-1)], at = 1:2,col=c("3","4"),boxwex = 0.6,xlim=c(0,9),xaxt = "n")
  abline(h=gamma,lty = 2)
  boxplot(pm[,4:5+2*(k-1)], at = 1:2+3 ,col=c("3","4"),add=T,boxwex = 0.6,xaxt = "n")
  boxplot(pm[,6:7+2*(k-1)], at = 1:2+6,col=c("3","4"),add=T,boxwex = 0.6,xaxt = "n")
  axis(1, at = seq(1.5,9,3), labels = c("MAF=0.05","MAF=0.10","MAF=0.15"), tick = TRUE)
  axis(1, at = seq(1.5,9,3), labels = c("MAF=0.20","MAF=0.25","MAF=0.30"), tick = TRUE)
  axis(1, at = seq(1.5,9,3), labels = c("MAF=0.35","MAF=0.40","MAF=0.45"), tick = TRUE)
}



###maf
# setwd("C:/Users/o0/Desktop/paper/0217/maf-g/")
setwd("/home/o0/Desktop/paper/0307/J=3/maf-x/maf-x2000/")
e=c("maf1.txt","maf2.txt","maf3.txt","maf4.txt","maf5.txt","maf6.txt","maf7.txt","maf8.txt","maf9.txt")
ma=seq(0.05,0.45,0.05)
title=sapply(ma,function(x){paste("MAF=",x)})
beta=log(c(1.1,1.2,1.3,1.4,1.5))
# beta=sort(c(-beta,0,beta))
par(mfrow=c(1,1))
type1=NULL
# par(mfrow=c(1,1))
for (i in 1:9){
  result=read.table(e[i],header = T,sep = ",")[,-1]
  type1=rbind(type1,result[1,])
  result=result[2:6,7:9]
  ###power
  plot(beta,result[,1],type="o",pch=1,xlim = c(log(1.1),log(1.55)),ylim = c(0,1),col=3,xlab = expression(beta),ylab="Power",main=(title[i]),xaxt="n")
  axis(1,at=c(log(1.1),log(1.2),log(1.3),log(1.4),log(1.5)),c("ln1.1","ln1.2","ln1.3","ln1.4","ln1.5"))
  lines(beta,result[,2],type="o",col=4)
  lines(beta,result[,3],type="o",col=2,pch=6)
  
  legend(0.36,0.25, c("pT","mT","wT"), lty = 1,pch=c(1,1,6),col=c(3,4,2),cex=0.8,text.font = 0.2)

}

##beta


# re=read.table("result",header = T)
write.table(cbind(ma,type1[,7:9]), file = "f.txt", sep = "&", row.names=FALSE, col.names=FALSE,quote = F)



result=read.table(e[i],header = T,sep = ",")[,-1]
result=result[-1,]
#####sqrt of mse plot
plot(beta,sqrt(result[,4]),type="o",xlim = c(-log(1.55),log(1.55)),ylim=c(0.085,0.12),col=3,xlab = expression(beta),ylab="The Square of MSE",xaxt="n",main=(title[i]))
axis(1,at=c(log(1.1),log(1.2),log(1.3),log(1.4),log(1.5)),c("ln1.1","ln1.2","ln1.3","ln1.4","ln1.5"))
axis(1,at=beta,c("-ln1.5","-ln1.4","-ln1.3","-ln1.2","-ln1.1","0","ln1.1","ln1.2","ln1.3","ln1.4","ln1.5"))
lines(beta,sqrt(result[,5]),type="o",col=4)
lines(beta,sqrt(result[,6]),type="o",col=2)
legend(0.1,0.08, c("pro","mod","new"), col=c(3,4,2),lty = 1,pch=c(1,1,1),cex=0.8)
