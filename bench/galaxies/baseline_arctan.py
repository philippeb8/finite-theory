"""Empirical baseline with the same parameter count and objective as the survey fit: v = (2/pi) V arctan(x/rt) + omega x,
omega set so the curve passes through the last point; V and rt fitted by minimising the standard deviation of the residuals.
Usage: python3 baseline_arctan.py <grc folder>  (after galaxy_fit.py; adds std_arctan to galaxies.csv)"""
import numpy as np, json, sys, os, glob, csv
from scipy.optimize import minimize
GRC=sys.argv[1]; KPC=3.08567758128e19
vel={}
for l in open(os.path.join(GRC,'db/0/grc.txt')):
    f=l.split(); vel.setdefault(f[0].upper(),[]).append((float(f[2])*KPC,float(f[3])*1000))
for p in sorted(glob.glob(os.path.join(GRC,'db/1/**/rc.dat'),recursive=True)):
    n=os.path.basename(os.path.dirname(p)).upper()
    for l in open(p):
        f=l.split()
        if len(f)>=2: vel.setdefault(n,[]).append((float(f[0])*KPC,float(f[1])*1000))
R=json.load(open('galaxies.json'))
def sd(p,x,y):
    V,rt=np.exp(p); v=2/np.pi*V*np.arctan(x/rt); w=(y[-1]-v[-1])/x[-1]; return np.std(y-(v+w*x))
for r in R:
    x=np.array([a for a,b in vel[r['name']]]); y=np.array([b for a,b in vel[r['name']]])
    best=None
    for V0 in (0.5,1,2):
        for f in (0.03,0.1,0.3,1):
            res=minimize(sd,np.log([V0*max(y.max(),1),f*x.max()]),args=(x,y),method='Nelder-Mead',options=dict(xatol=1e-4,fatol=1e-2,maxiter=2000))
            if best is None or res.fun<best.fun: best=res
    r['std_arctan']=float(best.fun)
json.dump(R,open('galaxies.json','w'))
with open('galaxies.csv','w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(R[0])); wr.writeheader(); [wr.writerow(q) for q in R]
a=np.array([q['std_arctan'] for q in R]); s=np.array([q['std'] for q in R])
print('arctan baseline mean std %.0f m/s; FT %.0f; FT better in %d of %d'%(a.mean(),s.mean(),(s<a).sum(),len(R)))
