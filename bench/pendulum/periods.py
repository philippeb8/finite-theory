"""Periods of the pendulum runs (Section 5.2, Fig. 10 a) from the tracked bob positions.
Usage: python3 periods.py            (reads tracks/*.csv, writes periods.json)
Each track: time (s), x (px), y (px), blob size; produced by track_pendulum.py from the videos.
Both records are used over the same swings: the 21 turning points at the near end and the 20 at the far end that
follow the release. The period is the slope of the turning times against the swing count.
The amplitude of each swing is followed by the peak-to-peak height of the bob in the image, scaled to the release
amplitude; the period at zero amplitude, T0, removes the amplitude dependence 1 + th^2/16 + 11 th^4/3072 swing by swing.
The crossings of the middle of the swing give a second estimate of the period."""
import numpy as np, json, os
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
fps=60000/1001
here=os.path.dirname(os.path.abspath(__file__))
runs={'IMG_8754_steel':(2.0,60.4),'IMG_8755_magnets':(2.0,60.4)}   # swinging part of each record (s), the same for both
TH_REL=np.radians(20.49)                                           # released 70 cm from the vertical on the 2 m rope
amp=lambda th: 1+th**2/16+11*th**4/3072                            # period / period at zero amplitude
def slope(tp,c=None):
    c=np.arange(len(tp)) if c is None else c; f=np.polyfit(c,tp,1); r=tp-np.polyval(f,c)
    return f[0], r.std(ddof=2)/np.sqrt(((c-c.mean())**2).sum())
out={}
for n,(t0,t1) in runs.items():
    t,x,y,s=np.loadtxt(os.path.join(here,'tracks',n+'.csv'),delimiter=',').T
    ok=np.isfinite(y); y=np.interp(t,t[ok],y[ok]); yf=median_filter(y,size=9)
    sel=(t>=t0)&(t<=t1); tt,yy=t[sel],yf[sel]; res={}
    for sign,side in ((1,'near'),(-1,'far')):
        pk,_=find_peaks(sign*yy,distance=int(2.0*fps),prominence=60); tp=[]; ye=[]
        for p in pk:                                   # turning time and height refined with a parabola over +-0.25 s
            if p<15 or p+16>len(tt): continue
            w=slice(p-15,p+16); c=np.polyfit(tt[w]-tt[p],sign*yy[w],2)
            tp.append(tt[p]-c[1]/(2*c[0])); ye.append(sign*(c[2]-c[1]**2/(4*c[0])))
        res[side]=dict(turning_times=tp,turning_y=ye)
    tn,yn=np.array(res['near']['turning_times']),np.array(res['near']['turning_y'])
    tf,yfa=np.array(res['far']['turning_times']),np.array(res['far']['turning_y'])
    # peak-to-peak height of each half swing, its linear trend, and the amplitude scaled to the release
    pp=[];tm=[]
    for k in range(len(tf)):
        pp.append(yn[k]-yfa[k]); tm.append((tn[k]+tf[k])/2)
        if k+1<len(tn): pp.append(yn[k+1]-yfa[k]); tm.append((tn[k+1]+tf[k])/2)
    cpp=np.polyfit(tm,pp,1)
    for side in ('near','far'):
        tp=np.array(res[side]['turning_times']); T,e=slope(tp)
        t_rel=tn[0]-T/2                                # release: half a swing before the first near turning point
        th=lambda tq: TH_REL*np.polyval(cpp,tq)/np.polyval(cpp,t_rel)
        cum=np.concatenate([[0],np.cumsum(amp(th((tp[1:]+tp[:-1])/2)))])
        T0,e0=slope(tp,cum)
        res[side].update(period=T,error=e,period_zero_amplitude=T0,error_zero_amplitude=e0)
    # crossings of the middle of the swing, between consecutive turning points
    ev=sorted([(a,b,'n') for a,b in zip(tn,yn)]+[(a,b,'f') for a,b in zip(tf,yfa)])
    yc=median_filter(y,size=5); cr={'n':[],'f':[]}
    for (ta,ya,sa),(tb,yb,sb) in zip(ev[:-1],ev[1:]):
        m=(t>ta)&(t<tb); q=yc[m]-(ya+yb)/2; tq=t[m]; i=np.where(np.sign(q[:-1])!=np.sign(q[1:]))[0]
        if len(i)==0: continue
        i=i[len(i)//2]; w=slice(max(i-4,0),i+6); c=np.polyfit(tq[w],q[w],1); cr[sa].append(-c[1]/c[0])
    Tc=[slope(np.array(cr[k]))[0] for k in ('n','f')]
    T=np.mean([res[q]['period'] for q in ('near','far')]); e=np.hypot(res['near']['error'],res['far']['error'])/2
    T0=np.mean([res[q]['period_zero_amplitude'] for q in ('near','far')])
    e0=np.hypot(res['near']['error_zero_amplitude'],res['far']['error_zero_amplitude'])/2
    out[n]=dict(res,mean_period=T,mean_error=e,period_zero_amplitude=T0,error_zero_amplitude=e0,
                crossing_period=float(np.mean(Tc)),peak_to_peak=dict(t=tm,px=pp,trend=cpp.tolist()),
                swings=len(tn)-1,turning_points=len(tn)+len(tf),
                end_amplitude_deg=float(np.degrees(TH_REL*np.polyval(cpp,tn[-1])/np.polyval(cpp,tn[0]-T/2))))
    print('%-18s T = %.4f +- %.4f s over %d swings (near %.4f, far %.4f); %d turning points; crossings %.4f s;'
          ' amplitude 20.5 -> %.1f deg; zero amplitude %.4f s'%(n,T,e,out[n]['swings'],res['near']['period'],res['far']['period'],
          out[n]['turning_points'],out[n]['crossing_period'],out[n]['end_amplitude_deg'],T0))
S,M=out['IMG_8754_steel'],out['IMG_8755_magnets']
d=M['mean_period']/S['mean_period']-1; ed=np.hypot(S['mean_error']/S['mean_period'],M['mean_error']/M['mean_period'])
d0=M['period_zero_amplitude']/S['period_zero_amplitude']-1
print('magnets longer by %.3f +- %.3f %% (crossings %.3f %%); at zero amplitude %.3f %%'%(100*d,100*ed,
      100*(M['crossing_period']/S['crossing_period']-1),100*d0))
json.dump(out,open(os.path.join(here,'periods.json'),'w'),indent=1)
