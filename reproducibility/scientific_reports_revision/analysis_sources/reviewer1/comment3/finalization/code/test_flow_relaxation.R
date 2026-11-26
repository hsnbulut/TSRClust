# Exhaustive small aggregate-count tables,not simulated observations or fitted models.
f<-grep('^--file=',commandArgs(),value=TRUE)[1];F<-dirname(dirname(normalizePath(sub('^--file=','',f))))
p<-parse(file.path(F,'reproduce_final_bounds.R'))
ix<-which(vapply(p,function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('flow_bounds')),logical(1)))
stopifnot(length(ix)==1);eval(p[[ix]])
comp<-function(n){a<-expand.grid(P=0:n,V=0:n,W=0:n,Z=0:n);as.matrix(a[rowSums(a)==n,,drop=FALSE])}
ncheck<-0L
for(M in 1:3)for(U in 0:2)for(N in 0:2){
 ma<-comp(M);ua<-comp(U);na<-comp(N)
 for(i in seq_len(nrow(ma)))for(j in seq_len(nrow(ua)))for(k in seq_len(nrow(na))){
  a<-rbind(ma[i,],ua[j,],na[k,]);v<-sum(a[,'V']);w<-sum(a[,'W']);pp<-sum(a[,'P'])
  if((v>0&&a[2,'V']==0)||(pp>0&&a[1,'P']==0))next
  for(direct in c(FALSE,TRUE)){
   b<-flow_bounds(M,U,N,a[2,'V'],sum(a[c(1,3),'V']),a[1,'P'],sum(a[c(1,3),c('V','W')]),a[3,'Z'],sum(a[,'Z']),as.integer(w>0),as.integer(v>0),if(direct)sum(a[1,c('V','W')]) else NA_integer_)
   stopifnot(a[1,'V']>=b['LV'],a[1,'V']<=b['UV'],sum(a[1,c('V','W')])>=b['LE'],sum(a[1,c('V','W')])<=b['UE'])
  }
  ncheck<-ncheck+1L
 }
}
# Fixed-pair algebra,not a claim about re-optimized full-partition matching.
for(u in 1:10)for(tp in 0:u)for(fp in 0:10)stopifnot(2*tp/(u+tp+fp+1)<=2*tp/(u+tp+fp))
write.csv(data.frame(test=c('Independent flow bounds contain every admissible small fixed role-count table','Adding nonminority counts never increases fixed-pair F1'),n_tables_or_count_cases=c(ncheck,sum(2:11)*11),status='PASS',scope='Deterministic integer/count algebra only;no data generation or model fitting'),file.path(F,'evidence/independent_math_regressions.csv'),row.names=FALSE)
cat('INDEPENDENT_FLOW_AND_FIXED_F1_REGRESSION_OK',ncheck,'count tables;with and without direct E\n')
