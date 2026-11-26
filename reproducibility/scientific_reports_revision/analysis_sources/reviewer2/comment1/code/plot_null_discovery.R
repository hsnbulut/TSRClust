f<-grep('^--file=',commandArgs(),value=TRUE)[1];C<-dirname(dirname(normalizePath(sub('^--file=','',f))))
x<-read.csv(file.path(C,'tables/false_micro_discovery_by_noise.csv'))
families<-c('uniform','correlated_gaussian','heterogeneous','heavy_tail','bridge','correlated_sheet')
labels<-c('Uniform','Corr. Gaussian','Beta gradient','Multiv. t3','Bridge','Corr. sheet')
pdf(file.path(C,'figures/false_discovery_by_noise.pdf'),width=10,height=7.7,pointsize=11)
layout(matrix(c(1,2,3,4,5,5),3,byrow=TRUE),heights=c(1,1,.25))
par(mar=c(4.1,7.6,2.3,1),las=1)
for(i in 1:4){
 p<-c(2,2,5,5)[i];alpha<-c(.1,.2,.1,.2)[i];lim<-if(p==2)50 else 2
 plot(c(-lim*.015,lim),c(.5,6.5),type='n',yaxt='n',xlab='Native false-discovery frequency (%)',ylab='',main=sprintf('p = %s, trimming = %.2f, residual n = %s',p,alpha,500*alpha),bty='l')
 axis(2,at=6:1,labels=labels);abline(v=pretty(c(0,lim)),col='#EEEEEE')
 for(j in 1:2){
  a<-x[x$experiment=='A-END'&x$p==p&x$requested_alpha==alpha&x$rule==c('legacy_rule','pvalue_only')[j],];a<-a[match(families,a$noise_type),];yy<-(6:1)+(j-1.5)*.16;co<-c('#235D8A','#B65C24')[j]
  arrows(100*a$native_lower,yy,100*a$native_upper,yy,angle=90,code=3,length=.025,col=co)
  points(100*a$native_rate,yy,col=co,pch=c(16,17)[j],cex=.8)
 }
}
par(mar=c(0,0,0,0));plot.new();legend(.5,.8,c('Legacy rule','P-value only'),col=c('#235D8A','#B65C24'),pch=c(16,17),bty='n',horiz=TRUE,xjust=.5,yjust=.5)
text(.5,.3,'1,000 independent datasets per fixed cell. Bars: pointwise exact binomial 95% intervals. Panel scales differ.',cex=.82)
text(.5,.05,'Native Stage-II acceptance is distinct from external eligible extras; correlated-sheet p=5, trim=.20 has external 40/1,000 but native zero.',cex=.7)
dev.off()
