# Valid aggregate-count bounds; no ARI/F1 inversion and no membership recovery.
origin_bounds <- function(M,U,N,TP,FP,FN,macroP,FPE,noise0,all0,spurious) {
 nonret<-M-macroP
 min_minor_unmatched<-pmax(0,spurious-(FPE-FP))
 max_minor0<-FN-min_minor_unmatched
 stopifnot(all(max_minor0>=0))
 LE<-pmax(0,FPE-(N-noise0),nonret+noise0-all0)
 UE<-pmin(nonret,FPE,nonret+noise0-all0+max_minor0)
 LV<-pmax(0,FP-(N-noise0),LE-(FPE-FP));UV<-pmin(FP,UE)
 LE<-pmax(LE,LV);UE<-pmin(UE,UV+(FPE-FP))
 LV<-pmax(LV,LE-(FPE-FP));UV<-pmin(UV,UE)
 stopifnot(all(LE<=UE),all(LV<=UV))
 list(LE=LE,UE=UE,LV=LV,UV=UV)
}
