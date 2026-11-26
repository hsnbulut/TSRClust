# Portable entry point for the unchanged submitted-summary arithmetic below.
# Historical raw-dataset/cross-installation checks are excluded, not claimed passed.
a<-commandArgs(TRUE);stopifnot(length(a)==2)
G<-normalizePath(a[1]);C<-normalizePath(a[2],mustWork=FALSE)
dir.create(C,recursive=TRUE,showWarnings=FALSE)
out<-file.path(C,'tables');ev<-file.path(C,'evidence')
dir.create(out,showWarnings=FALSE);dir.create(ev,showWarnings=FALSE)
write<-function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na='')
get<-function(name)read.csv(file.path(G,'results/real_data_results',name),stringsAsFactors=FALSE)
comp<-get('shuttle_microcluster_composition.csv')
cls<-subset(get('class_trimming_rates.csv'),dataset=='Statlog_Shuttle')
fit<-subset(get('dataset_fit_summary.csv'),dataset=='Statlog_Shuttle')
res<-subset(get('residual_external_eval.csv'),dataset=='Statlog_Shuttle'&method=='TSRClust_validated')
classes<-cls$local_class
original<-file.path(G,'reproducibility')
# Verify how these summaries were computed: unrounded integer tabulations of saved assignments.
script <- readLines(file.path(original,'real_data/02_realdata_v2.R'),warn=FALSE)
stopifnot(any(grepl('tt <- table(truth[idx])',script,fixed=TRUE)),
          any(grepl('purity=max(tt)/length(idx)',script,fixed=TRUE)),
          any(grepl('cluster_purity(local_truth[r0], final[r0], positive_only=TRUE)',script,fixed=TRUE)))
# At sizes <=114, one contaminant changes purity by >=1/114, far above CSV precision.
# Thus size + dominant class + purity=1 are sufficient counts, not recovered row memberships.
b <- which(cls$local_class=='Bypass'); total_b <- cls$total[b]; trim_b <- cls$trimmed[b]
comp$Bypass <- comp$size
for (cl in classes[classes!='Bypass']) comp[[cl]] <- 0L
comp$other_classes <- 0L
comp$accepted_share <- comp$size/sum(comp$size)
comp$original_Bypass_share <- comp$size/total_b
comp$trimmed_Bypass_share <- comp$size/trim_b
stopifnot(all(rowSums(comp[,classes])==comp$size),sum(comp$Bypass)==511)
write(comp,'shuttle_microcluster_composition.csv')

cls$StageI_macro <- cls$total-cls$trimmed
cls$accepted_StageII <- ifelse(cls$local_class=='Bypass',511L,0L)
cls$final_zero <- cls$trimmed-cls$accepted_StageII
cls$macro_fraction_original <- cls$StageI_macro/cls$total
cls$accepted_fraction_original <- cls$accepted_StageII/cls$total
cls$accepted_fraction_trimmed <- cls$accepted_StageII/cls$trimmed
cls$zero_fraction_original <- cls$final_zero/cls$total
cls$zero_fraction_trimmed <- cls$final_zero/cls$trimmed
stopifnot(all(cls$StageI_macro+cls$accepted_StageII+cls$final_zero==cls$total),
          all(cls$accepted_StageII+cls$final_zero==cls$trimmed),
          sum(cls$final_zero)==2389,cls$StageI_macro[b]==958,cls$final_zero[b]==1798)
write(cls,'shuttle_class_flow.csv')

tab <- matrix(0L,7,12,dimnames=list(classes,c('0',as.character(comp$cluster))))
tab[,'0'] <- cls$final_zero;tab['Bypass',-1] <- comp$size
stopifnot(all(rowSums(tab)==cls$trimmed),sum(tab)==2900,all(colSums(tab)[-1]==comp$size))
write(data.frame(external_class=rownames(tab),tab,check.names=FALSE),'shuttle_residual_contingency.csv')
choose2 <- function(n) n*(n-1)/2
T <- choose2(sum(tab)); A <- sum(choose2(tab)); R <- sum(choose2(rowSums(tab)))
P <- sum(choose2(colSums(tab))); E <- R*P/T; denominator <- (R+P)/2-E
ari <- (A-E)/denominator
stopifnot(abs(ari-res$ResidualLocalARI)<1e-14,abs(sum(comp$size)/2900-res$ResidualAssignedFraction)<1e-14)
accepted_split <- choose2(511)-sum(choose2(comp$size))
accepted_zero <- 511*1798
pair <- data.frame(quantity=c('All unordered residual pairs','Same external class and same predicted label',
 'Same external class pairs','Same predicted label pairs','Same class different labels',
 'Same label different classes','Different class and different labels',
 'Bypass pairs across accepted microclusters','Bypass pairs accepted versus final zero',
 'Chance expected same-both pairs','ARI denominator','ResidualLocalARI'),
 value=c(T,A,R,P,R-A,P-A,T-R-P+A,accepted_split,accepted_zero,E,denominator,ari))
stopifnot(accepted_split+accepted_zero==R-A,tab['Bypass','0']==1798)
write(pair,'shuttle_pair_decomposition.csv')

coverage <- data.frame(quantity=c('Original Bypass coverage','Trimmed Bypass coverage',
 'Overall residual assignment','Bypass retained in Stage-I macros','Trimmed Bypass still zero',
 'Original Bypass still zero','Non-Bypass residual share','Non-Bypass residual accepted'),
 numerator=c(511,511,511,958,1798,1798,591,0),
 denominator=c(3267,2309,2900,3267,2309,3267,2900,591))
coverage$fraction <- coverage$numerator/coverage$denominator
write(coverage,'shuttle_coverage_summary.csv')
write(data.frame(state=c('Stage-I macro assigned','Accepted Stage-II microclusters','Unassigned residual'),
                 count=c(958,511,1798),original_Bypass_total=3267,
                 fraction=c(958,511,1798)/3267),'shuttle_Bypass_distribution.csv')

cat('PUBLIC_SHUTTLE_COUNTS_OK; Bypass 511/3267 and 511/2309; ARI=',sprintf('%.16f',ari),'\n')
