"""Periods of the pendulum runs (Section 5.2, Fig. 10 a) from the tracked bob positions.
Usage: python3 periods.py            (reads tracks/*.csv, writes periods.json)
Each track: time (s), x (px), y (px), blob size; produced by track_pendulum.py from the videos."""
import numpy as np, json, os
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
fps=60000/1001
here=os.path.dirname(os.path.abspath(__file__))
runs={'IMG_8754_steel':(2.0,53.5),'IMG_8755_magnets':(2.0,60.3)}   # swinging part of each record (s)
out={}
for n,(t0,t1) in runs.items():
    t,x,y,s=np.loadtxt(os.path.join(here,'tracks',n+'.csv'),delimiter=',').T
    ok=np.isfinite(y); y=np.interp(t,t[ok],y[ok]); yf=median_filter(y,size=9)
    sel=(t>=t0)&(t<=t1); tt,yy=t[sel],yf[sel]; res={}
    for sign,side in ((1,'near'),(-1,'far')):
        pk,_=find_peaks(sign*yy,distance=int(2.0*fps),prominence=60); tp=[]
        for p in pk:                                   # turning time refined with a parabola over +-0.25 s
            w=slice(max(p-15,0),min(p+16,len(tt))); c=np.polyfit(tt[w]-tt[p],sign*yy[w],2); tp.append(tt[p]-c[1]/(2*c[0]))
        tp=np.array(tp); k=np.arange(len(tp)); f=np.polyfit(k,tp,1); r=tp-np.polyval(f,k)
        res[side]=dict(turning_times=tp.tolist(),period=f[0],error=r.std()/np.sqrt(((k-k.mean())**2).sum()))
    T=np.mean([res[q]['period'] for q in res]); e=np.sqrt(sum(res[q]['error']**2 for q in res))/2
    out[n]=dict(res,mean_period=T,mean_error=e,swings=len(res['near']['turning_times'])-1)
    print('%-18s T = %.4f +- %.4f s over %d swings (near %.4f, far %.4f)'%(n,T,e,out[n]['swings'],res['near']['period'],res['far']['period']))
Ts,Tm=out['IMG_8754_steel']['mean_period'],out['IMG_8755_magnets']['mean_period']
es,em=out['IMG_8754_steel']['mean_error'],out['IMG_8755_magnets']['mean_error']
d=Tm/Ts-1; ed=np.hypot(es/Ts,em/Tm)
print('magnets longer by %.3f +- %.3f %%'%(100*d,100*ed))
json.dump(out,open(os.path.join(here,'periods.json'),'w'),indent=1)
