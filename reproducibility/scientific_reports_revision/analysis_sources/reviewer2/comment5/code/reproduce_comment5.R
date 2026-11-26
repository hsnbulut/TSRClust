# Read-only arithmetic from submitted sufficient summaries. No estimators are called.
args <- commandArgs(TRUE)
project <- normalizePath(if (length(args)) args[1] else '.', mustWork=TRUE)
C <- file.path(project,'Revision1/reviewer2/comment5')
original <- file.path(project,'Scientific Reports Manuscript/analysis')
out <- file.path(C,'tables'); ev <- file.path(C,'evidence')
write <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,na='')
get <- function(name) read.csv(file.path(original,'v2/results/realdata',name),stringsAsFactors=FALSE)
comp <- get('shuttle_microcluster_composition.csv')
cls <- subset(get('class_trimming_rates.csv'),dataset=='Statlog_Shuttle')
fit <- subset(get('dataset_fit_summary.csv'),dataset=='Statlog_Shuttle')
res <- subset(get('residual_external_eval.csv'),dataset=='Statlog_Shuttle' & method=='TSRClust_validated')
tim <- subset(get('timings.csv'),dataset=='Statlog_Shuttle' & method=='TSRClust_validated')
global <- subset(get('external_macro_ari.csv'),dataset=='Statlog_Shuttle' & method=='TSRClust_validated')

# Load only the existing data object; do not load the package DLL or install anything.
datafile <- Sys.getenv('TSR_SHUTTLE_RDA')
e <- new.env(); load(datafile,envir=e); d <- e$Shuttle
classes <- levels(d[[10]])
stopifnot(nrow(d)==58000,ncol(d)==10,length(classes)==7,classes[5]=='Bypass',!anyNA(d),
          identical(classes,cls$local_class),all(as.integer(table(d[[10]]))==cls$total),
          sum(cls$total)==58000,sum(cls$trimmed)==2900,nrow(fit)==1,
          fit$n==58000,fit$p==9,fit$k==2,fit$alpha==.05,fit$c==20,
          fit$residual_n==2900,fit$micro_q==11,fit$micro_p_value==.02,
          identical(comp$cluster,3:13),sum(comp$size)==511,
          all(comp$dominant_class=='Bypass'),all(comp$purity==1),all(comp$dataset=='Statlog_Shuttle'))

# Verify how these summaries were computed: unrounded integer tabulations of saved assignments.
script <- readLines(file.path(original,'v2/scripts/02_realdata_v2.R'),warn=FALSE)
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

# Confirm the distinct historical fit cannot be substituted for this submitted result.
old <- readRDS(file.path(project,'Eski Dosyalar/Shuttle_TSRClust_results/Shuttle_TSRClust_analysis_summary.rds'))
stopifnot(old$selected_c==1,old$n_micro_clusters==3,identical(old$truth,d[[10]]))
write(data.frame(source=c('Submitted V2 result','Excluded historical RDS'),K=2,
                 c=c(20,old$selected_c),microclusters=c(11,old$n_micro_clusters),
                 residual_n=c(2900,old$n_r0),used_for_memberships=FALSE), 'shuttle_source_versions.csv')
orderfile <- file.path(ev,'qa/ordered_input.csv')
write.csv(d,orderfile,row.names=FALSE)
other <- new.env();load('${TSR_SHUTTLE_RDA}',envir=other)
stopifnot(identical(d,other$Shuttle))
write(data.frame(code=seq_along(classes),class=classes,count=as.integer(table(d[[10]]))), 'shuttle_data_classes.csv')

geometry <- data.frame(diagnostic=c('Selected epsilon','Candidate epsilon grid','Density/novelty scores',
 'Centroids/distances','Within-cluster dispersion','Observation-level raw/final memberships',
 'Number of selected epsilon values','Support','Reported global MC p value'),
 available=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE),
 value=c('Not exported','Not exported','Not exported','Not exported','Not exported','Not exported',
         'One per final solution by executed source','18 = max(5,2*9)','0.02'),
 basis=c(rep('Submitted output inventory',6),'core_tsr.R returns only best$labels and best$eps',
         'core_tsr.R and p=9','dataset_fit_summary.csv'))
write(geometry,'shuttle_geometry_availability.csv')
write(data.frame(metric=c('ExternalLabelARI','ResidualLocalARI','ResidualAssignedFraction',
 'preprocess_seconds','macro_seconds','micro_seconds','end_to_end_seconds'),
 value=c(global$MacroLabelARI,ari,511/2900,tim$preprocess_seconds,tim$macro_seconds,
         tim$micro_seconds,tim$end_to_end_seconds)), 'shuttle_original_metrics.csv')
checks <- data.frame(check=c('Full dataset and class mapping','Original fit settings','Composition integer sums',
 'Class-flow conservation','Exact residual contingency','Archived ARI reproduction',
 'Residual assignment reproduction','Historical source exclusion','Two local data versions identical'),
 status='PASS',detail=c('58000 rows,9 predictors,7 levels;Bypass code5;all source totals match',
 'K2 alpha.05 c20,2900 residuals,11 clusters,p=.02','11 sizes sum511;all7 class counts sum to size',
 'Every class:macro+accepted+zero=original','7x12 counts sum2900;zero column2389',
 sprintf('Computed %.16f;stored %.16f',ari,res$ResidualLocalARI),'511/2900 matches stored fraction',
 'Old c1,3-cluster fit rejected;truth row order matches','R4.2 andR4.5 Shuttle objects identical'))
write.csv(checks,file.path(ev,'analysis_checks.csv'),row.names=FALSE)
pdf(file.path(C,'figures/shuttle_microcluster_sizes.pdf'),width=7,height=5.3,useDingbats=FALSE)
par(mar=c(6,4,3,1),family='Helvetica')
bar <- barplot(comp$size,names.arg=comp$cluster,col='#466780',border=NA,
               ylim=c(0,130),xlab='',ylab='Observations',
               main='Eleven accepted Shuttle microclusters',las=1)
text(bar,comp$size+4,comp$size,cex=.85)
mtext('Native accepted microcluster label',side=1,line=2.5,cex=.9)
mtext('Submitted summaries: all 511 observations are Bypass; no refit',side=1,line=4.5,cex=.75)
dev.off()
pdf(file.path(C,'figures/shuttle_class_coverage.pdf'),width=8,height=5.3,useDingbats=FALSE)
par(mar=c(6,12,3,2),family='Helvetica')
counts <- c(1798,511,958); labels <- c('Unassigned residual','Accepted microclusters','Stage-I macro assigned')
bar <- barplot(counts,names.arg=labels,horiz=TRUE,col=c('#999999','#466780','#A7B6C1'),
               border=NA,xlim=c(0,2600),xlab='',
               main='Where the original Bypass class remains',las=1)
text(counts+35,bar,sprintf('%s (%.2f%%)',format(counts,big.mark=','),100*counts/3267),pos=4,cex=.85)
mtext('Original Bypass observations (denominator 3,267)',side=1,line=2.5,cex=.85)
mtext('Stage-I trimmed Bypass: 2,309; accepted 511/2,309 = 22.13%',side=1,line=4.5,cex=.75)
dev.off()
cat('R2C5_READ_ONLY_AUDIT_OK; sizes=',paste(comp$size,collapse=','),'; ARI=',sprintf('%.16f',ari),'; NO_FITS_NO_JOBS\n')
