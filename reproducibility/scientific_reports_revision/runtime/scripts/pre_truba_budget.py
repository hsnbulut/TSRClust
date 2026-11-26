"""Pilot-based planning projections; not hardware-confirmed TRUBA predictions."""
from pathlib import Path
import csv,math,statistics
r=Path(__file__).resolve().parents[1]
def read(path):return list(csv.DictReader((r/path).open()))
design=read('config/design.csv');raw=read('results/pilot_R452_pre_truba/replications.csv');pilot=[x for x in raw if x['status']=='OK']
def kind(z):return 'selected' if z['regime']=='data_selected' else 'fixed'
rows=[]
for d in design:
 pool=[x for x in pilot if x['method']==d['method'] and x['mode']==d['mode'] and kind(x)==kind(d) and (d['method']!='TSRClust' or x['reference']==d['reference'])]
 exact=[x for x in pool if x['cell_id']==d['cell_id']]
 if exact:sec=statistics.mean(float(x['end_to_end_seconds']) for x in exact);basis='same cell, two pilot replications'
 elif pool:
  distances=[abs(math.log(float(x['n'])/float(d['n'])))+abs(math.log(float(x['p'])/float(d['p']))) for x in pool];md=min(distances)
  near=[x for x,dist in zip(pool,distances) if dist<=md+1e-9];values=[]
  for x in near:
   scale=max(1,float(d['n'])/float(x['n']))*max(1,float(d['p'])/float(x['p']))
   seconds=float(x['end_to_end_seconds'])
   if d['method']=='TSRClust':seconds+=float(x['micro_seconds'])*(max(1,int(d['B'])/int(x['B']))-1)
   if kind(d)=='selected':scale*=max(1,int(d['kmax'])/int(x['kmax']))
   values.append(seconds*scale)
  sec=max(values);basis='maximum nearest n,p pilot; B micro-time and K-range scaled proxy'
 else:sec=None;basis='unmeasured'
 rows.append(dict(cell_id=d['cell_id'],tier=d['tier'],experiment=d['experiment'],method=d['method'],reps=d['reps'],seconds_per_fit='' if sec is None else sec,CPU_hours='' if sec is None else sec*int(d['reps'])/3600,basis=basis))
with (r/'results/pre_truba/workload_pilot_projection.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
summary=[]
for tier in ['1','2']:
 z=[x for x in rows if x['tier']==tier];sec=sum(float(x['CPU_hours']) for x in z if x['CPU_hours']!='')
 summary.append(dict(tier=tier,cells=len(z),fits=sum(int(x['reps']) for x in z),proxy_CPU_hours=sec,unmeasured_fits=sum(int(x['reps']) for x in z if x['CPU_hours']==''),reserve_2x_CPU_hours=sec*2,reserve_4x_CPU_hours=sec*4))
with (r/'results/pre_truba/workload_budget_summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=summary[0].keys());w.writeheader();w.writerows(summary)
for row in summary:print(row)
