"""Fits of the floating-magnet drifts (Section 5.3, Table 5): constant force from a free onset time,
opposed by a water drag linear in the speed; steady speed = slope of the later part of the track, before the wall.
Usage: python3 fit_drift.py"""
import numpy as np, os
from scipy.optimize import least_squares
here=os.path.dirname(os.path.abspath(__file__))
# run: (track file, pixels per cm, sign of the drift direction in pixels, direction name)
runs=[('IMG_8733_afternoon_track.csv',880/57.0,-1,'north'),('IMG_8735_afternoon_track.csv',880/57.0,-1,'north'),
      ('IMG_8745_Double-South_track.csv',248/15.0,+1,'south'),('IMG_8746_Double-North_track.csv',248/15.0,+1,'south'),
      ('IMG_8747_Single-South_track.csv',248/15.0,+1,'south'),('IMG_8748_Single-North_track.csv',248/15.0,+1,'south')]
def model(p,t): x0,t0,vT,tau=p; s=np.clip(t-t0,0,None); return x0+vT*(s-tau*(1-np.exp(-s/tau)))
print('%-34s %-6s %10s %12s %10s %8s'%('run','toward','a0 (mm/s2)','steady (mm/s)','travel(cm)','onset(s)'))
for f,px,sg,dirn in runs:
    A=np.loadtxt(os.path.join(here,'tracks',f),delimiter=',')
    if 'afternoon' in f: t=A[:,0]; c=(A[:,1]+A[:,2])/2; w=A[:,2]-A[:,1]; tol=25     # time, left edge, right edge
    else: t,c,w=A[:,0],A[:,1],A[:,2]; tol=15                                        # time, centre, width
    ok=np.isfinite(c)&(np.abs(w-np.nanmedian(w))<tol); t,c=t[ok],c[ok]
    x=sg*(c-c[:10].mean())
    end=x.max(); sel=(x<end-15)&(t>t.max()*0.45); v=np.polyfit(t[sel],x[sel],1)[0]
    fs=x<end-15; best=None
    for t0 in (0,5,10,15,20):
        r=least_squares(lambda p: model(p,t[fs])-x[fs],[0,t0,v,5.0],bounds=([-50,-5,0,0.05],[50,60,50,200]))
        if best is None or r.cost<best.cost: best=r
    x0,t0,vT,tau=best.x
    print('%-34s %-6s %10.2f %12.1f %10.1f %8.1f'%(f.replace('_track.csv',''),dirn,vT/tau/px*10,v/px*10,end/px,t0))
