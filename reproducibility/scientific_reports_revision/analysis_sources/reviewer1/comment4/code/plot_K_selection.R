f<-grep('^--file=',commandArgs(),value=TRUE)[1];c4<-dirname(dirname(normalizePath(sub('^--file=','',f))));x<-read.csv(file.path(c4,'tables/performance_summary.csv'))
pdf(file.path(c4,'figures/K_selection_tradeoffs.pdf'),width=8,height=7,pointsize=10)
layout(matrix(c(1,2,3,4,5,5),nrow=3,byrow=TRUE),heights=c(1,1,.24));par(mar=c(4.1,4.2,2,1),las=1)
cols<-c('#22577A','#A34827');scs<-c('negative','positive','two_minorities','weak');titles<-c('Negative control','Isolated minority','Two minorities','Weak minority')
for(i in seq_along(scs)){
 plot(c(3,7),c(.4,1),type='n',xaxt='n',xlab='Upper candidate count m',ylab='Mean MajorARI',main=titles[i],bty='l');axis(1,at=c(3,5,7))
 b<-x[x$scenario==scs[i]&x$experiment=='D1-E'&x$method=='TSRClust'&x$metric=='MajorARI',];abline(h=b$mean,col='#666666',lty=2)
 for(j in 1:2){a<-x[x$scenario==scs[i]&x$experiment=='D2-E'&x$method==c('TSRClust','tEIGEN_BIC')[j]&x$metric=='MajorARI',];a<-a[order(a$kmax),];lines(a$kmax,a$mean,type='b',col=cols[j],pch=c(16,17)[j]);arrows(a$kmax,a$lower,a$kmax,a$upper,angle=90,code=3,length=.035,col=cols[j])}
}
par(mar=c(0,0,0,0));plot.new();legend(.5,.83,legend=c('Selected TSRClust','BIC tEIGEN','Known-K TSRClust'),col=c(cols,'#666666'),lty=c(1,1,2),pch=c(16,17,NA),bty='n',horiz=TRUE,xjust=.5,yjust=.5,cex=.9)
text(.5,.25,'Each scenario uses the same 100 datasets across methods and ranges. Bars are pointwise Monte Carlo 95% intervals.',cex=.75)
dev.off()
