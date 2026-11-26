# Exhaustive fixed small role-count tables, not generated observations or fits.
f<-grep('^--file=',commandArgs(),value=TRUE)[1];code<-dirname(normalizePath(sub('^--file=','',f)))
source(file.path(code,'origin_bounds.R'))
compositions<-function(n){z<-expand.grid(P=0:n,V=0:n,W=0:n,Z=0:n);as.matrix(z[rowSums(z)==n,])}
ma<-compositions(3);mi<-compositions(2);no<-compositions(2);checked<-0L
for(i in seq_len(nrow(ma)))for(j in seq_len(nrow(mi)))for(k in seq_len(nrow(no))){
 a<-rbind(ma[i,],mi[j,],no[k,]);sp<-as.integer(sum(a[,'W'])>0)
 # V is an empty union or a positively matched union with true-minority overlap.
 if(sum(a[,'V'])>0&&a[2,'V']==0)next
 b<-origin_bounds(3,2,2,a[2,'V'],sum(a[c(1,3),'V']),2-a[2,'V'],a[1,'P'],sum(a[c(1,3),c('V','W')]),a[3,'Z'],sum(a[,'Z']),sp)
 stopifnot(a[1,'V']>=b$LV,a[1,'V']<=b$UV,sum(a[1,c('V','W')])>=b$LE,sum(a[1,c('V','W')])<=b$UE)
 checked<-checked+1L
}
write.csv(data.frame(test='Exhaustive small role-count tables contain their known macro V/E counts',n_tables=checked,status='PASS'),file.path(dirname(code),'evidence/origin_bound_regression.csv'),row.names=FALSE)
cat('ORIGIN_BOUND_REGRESSION_OK',checked,'fixed count tables; no fitting\n')
