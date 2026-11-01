"""Summary of the survey fit (Sec. 4.3.1). Usage: python3 summary.py   (after galaxy_fit.py and baseline_arctan.py)"""
import numpy as np, csv, re, collections
R=list(csv.DictReader(open('galaxies.csv'))); f=lambda k: np.array([float(r[k]) for r in R])
sd,rms,sdN,a,h,p=f('std'),f('rms'),f('std_newton'),f('std_arctan'),f('h'),f('points'); src=[r['source'] for r in R]
def norm(n):
    for pre,rx in (('UGC',r'^(UGC|U|Z)0*(\d+)$'),('NGC',r'^(NGC|N)0*(\d+)$'),('IC',r'^(IC)0*(\d+)$')):
        m=re.match(rx,n)
        if m: return pre+m.group(2)
    return n
g=collections.defaultdict(list); [g[norm(r['name'])].append(r['name']) for r in R]; d={k:v for k,v in g.items() if len(v)>1}
print('curves %d (SPARC %d, Sofue atlas %d); objects appearing twice %d; galaxies %d'%(len(R),src.count('SPARC'),src.count('Sofue'),len(d),len(g)))
print('duplicates:',', '.join('/'.join(v) for v in d.values()))
print('points per curve: median %d, range %d to %d'%(np.median(p),p.min(),p.max()))
print('Finite Theory, Eq. (91): mean std of residuals %.0f m/s (median %.0f); mean rms %.0f m/s'%(sd.mean(),np.median(sd),rms.mean()))
print('empirical arctangent curve, same objective and parameter count: mean std %.0f m/s (median %.0f)'%(a.mean(),np.median(a)))
print('Newtonian limit of Eq. (91), lambda_g -> infinity: mean std %.0f m/s'%sdN.mean())
print('Finite Theory better than the arctangent curve on %d of %d curves'%((sd<a).sum(),len(R)))
print('fitted lambda_g: quartiles %.2e %.2e %.2e kg/m; at the grid bounds %d'%(*np.percentile(h,[25,50,75]),((h<=1.0000001e18)|(h>=9.9e21)).sum()))
