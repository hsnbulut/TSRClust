f<-grep('^--file=',commandArgs(),value=TRUE)[1];c3<-dirname(dirname(normalizePath(sub('^--file=','',f))));x<-read.csv(file.path(c3,'tables/scenario_metric_summaries.csv'))
pdf(file.path(c3,'figures/otrimle_macro_minority_tradeoff.pdf'),width=8,height=5,pointsize=10)
layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),heights=c(1,.18))
par(mar=c(4.4,4.2,2,1),las=1)
methods<-c('OTRIMLE_fixedK','TSRClust_validated','TCLUST');cols<-c('#A34827','#22577A','#555555')
for(cohort in c('Original architecture M1-M4','Original stress')){
 plot(c(0,1),c(0,1),type='n',xlab='Mean MajorARI',ylab='Mean MinorityClusterF1',main=if(grepl('architecture',cohort))'Original architecture' else 'Original stress',bty='l')
 rect(0,.7,.1,1,col='#F2F2F2',border=NA)
 for(j in seq_along(methods)){
  a<-x[x$cohort==cohort&x$method==methods[j]&x$metric=='MajorARI',c('cell_id','mean')]
  b<-x[x$cohort==cohort&x$method==methods[j]&x$metric=='MinorityClusterF1',c('cell_id','mean')]
  z<-merge(a,b,by='cell_id');points(z$mean.x,z$mean.y,col=cols[j],pch=c(16,4,1)[j],cex=.75)
 }
}
par(mar=c(0,0,0,0));plot.new()
legend(.5,.82,legend=c('OTRIMLE','TSRClust','TCLUST'),col=cols,pch=c(16,4,1),bty='n',cex=.85,horiz=TRUE,xjust=.5,yjust=.5)
text(.5,.28,'Each point is an existing design-cell mean. Grey region: descriptive low-macro/high-minority flag, not a selected validation sample.',cex=.67)
dev.off()
