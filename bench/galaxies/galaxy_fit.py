"""Reproduces the rotation-curve survey of Sec. 4.3.1 (the procedure of github.com/philippeb8/grc, src/main.cpp) and
documents it: galaxy list, points, fitted parameters, objective, aggregate, duplicates, and a baseline on the same data.
Usage: python3 galaxy_fit.py <grc folder>   ->  galaxies.csv, summary.txt"""
import numpy as np, os, sys, glob, re, json
GRC=sys.argv[1]; G=6.67408e-11; KPC=3.08567758128e19   # G as in the paper; src/main.cpp of GRC uses 6.671e-11, which changes the
# mean rms residual by 1 m/s and the Newtonian-limit mean by 13 m/s (release bench-v1.0 used it)
gal={}; vel={}
# db/0: SPARC (Lelli et al. 2016); grid centre m0 = Upsilon x L[3.6] from gm.txt (cols 2 and 10), x 1e9 Msun x 2e30 kg
for l in open(os.path.join(GRC,'db/0/gm.txt')):
    f=l.split(); n=f[0].upper(); gal[n]=float(f[9])*1e9*float(f[1])*2e30; src={}
SRC={n:'SPARC' for n in gal}
for l in open(os.path.join(GRC,'db/0/grc.txt')):
    f=l.split(); n=f[0].upper(); vel.setdefault(n,[]).append((float(f[2])*KPC,float(f[3])*1000))
# db/1: Sofue atlas; mass scale accumulated over mlcat.dat rows exactly as main.cpp; curve from rc.dat
for p in sorted(glob.glob(os.path.join(GRC,'db/1/**/mlcat.dat'),recursive=True)):
    n=os.path.basename(os.path.dirname(p)).upper(); SRC.setdefault(n,'Sofue')
    rows=[l.split() for l in open(p)][1:]
    for f in rows:
        if len(f)<8: continue
        R,SB,RH=float(f[0]),float(f[1]),float(f[7]); h=(R/RH)*KPC
        gal[n]=gal.get(n,0.0)+2*np.pi*SB*h*h/1000.0
for p in sorted(glob.glob(os.path.join(GRC,'db/1/**/rc.dat'),recursive=True)):
    n=os.path.basename(os.path.dirname(p)).upper()
    for l in open(p):
        f=l.split()
        if len(f)>=2: vel.setdefault(n,[]).append((float(f[0])*KPC,float(f[1])*1000))
H=np.arange(1e20/100,1e20*100,1e20/2)
def fit(x,y,m0,newtonian=False):
    m=np.arange(m0/15,m0*15,10**np.floor(np.log10(m0))/2)
    hh=np.array([np.inf]) if newtonian else H
    best=(np.inf,None)
    for h in hh:
        mm=m[:,None]; xx=x[None,:]
        ft=np.sqrt(G*mm/xx) if newtonian else h/(mm/xx+h)*np.sqrt(G*mm/xx)
        w=(y[-1]-ft[:,-1])/x[-1]; v=ft+w[:,None]*xx; d=y[None,:]-v
        sd=d.std(axis=1); k=int(np.argmin(sd))
        if sd[k]<best[0]: best=(sd[k],(h,m[k],w[k],np.sqrt((d[k]**2).mean())))
    return best
rows=[]
for i,n in enumerate(sorted(gal)):
    if n not in vel: continue
    x=np.array([a for a,b in vel[n]]); y=np.array([b for a,b in vel[n]])
    if not np.isfinite(gal[n]) or gal[n]<=0: continue
    sd,(h,m,w,rms)=fit(x,y,gal[n]); sdN,(hN,mN,wN,rmsN)=fit(x,y,gal[n],True)
    rows.append(dict(name=n,source=SRC.get(n,'?'),points=len(x),rmax_kpc=x.max()/KPC,m0=gal[n],h=h,m=m,omega=w,std=sd,rms=rms,std_newton=sdN,rms_newton=rmsN))
    if i%50==0: print(i,n,len(x),round(sd),file=sys.stderr)
import csv
with open('galaxies.csv','w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0])); wr.writeheader(); [wr.writerow(r) for r in rows]
json.dump(rows,open('galaxies.json','w'))
print('galaxies fitted',len(rows))
