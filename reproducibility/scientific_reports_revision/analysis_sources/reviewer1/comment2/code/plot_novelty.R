f<-grep('^--file=',commandArgs(),value=TRUE)[1];c2<-dirname(dirname(normalizePath(sub('^--file=','',f))))
x<-read.csv(file.path(c2,'tables/proximity_summary.csv'))
pdf(file.path(c2,'figures/novelty_proximity.pdf'),width=9,height=3.7,pointsize=10)
par(mfrow=c(1,3),mar=c(5.3,4.3,2.4,1),oma=c(1.4,0,0,0),las=1)
vs<-sort(unique(x$novelty));cols<-c('#22577A','#A34827')
for(prox in c('isolated','adjacent','overlapping')) {
 plot(1:5,rep(0,5),type='n',ylim=c(-.025,1),xaxt='n',xlab='',ylab=if(prox=='isolated')'Mean MinorityClusterF1' else '',main=prox,bty='l')
 axis(1,at=1:5,labels=sprintf('%.3f',vs),las=2,cex.axis=.85);abline(v=which(vs==.995),lty=3,col='gray50')
 for(j in 1:2){z<-x[x$proximity==prox&x$minority_ratio==c(.025,.05)[j],];z<-z[order(z$novelty),];lines(1:5,z$cluster_F1,col=cols[j],lty=j,pch=c(16,4)[j],type='b',lwd=1.4)}
 mtext('novelty probability',side=1,line=4)
 if(prox=='isolated')legend('bottomleft',legend=c('Minority fraction 0.025','Minority fraction 0.050'),col=cols,lty=1:2,pch=c(16,4),bty='n',cex=.8)
}
mtext('Existing B cells, legacy empirical reference and legacy rule; p = 5; 200 fits per point. Dotted line: default 0.995.',outer=TRUE,side=1,line=.1,cex=.8)
dev.off()
