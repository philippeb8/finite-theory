"""Fit of Eq. (91) to the rotation curve of the Milky Way (Secs. 2.2.1 and 4.3.1), by the procedure of galaxy_fit.py:
lambda_g on the same grid (1e18 to 1e22 kg/m), m between 1/15 and 15 times a mass scale, omega set so that the curve
passes through the last point, objective = standard deviation of the residuals about their mean (unweighted).
Data: milky_way_sofue2020.csv, the curve of Sofue (2020), 0.1 to 95.6 kpc. Mass scale: the stellar mass of the
Galaxy, 5e10 Msun (Bland-Hawthorn & Gerhard 2016, Sec. 2.2), with 2e30 kg per Msun as in galaxy_fit.py.
The Milky Way is not one of the 497 survey curves.
Principal fit: outside the bulge, from 2.5 kpc (the minimum of the curve, at 2.3 to 2.6 kpc). This inner limit was
chosen after inspecting the fits: over the whole curve the single central mass of Eq. (91) cannot follow the peak of
the bulge at 0.4 kpc. Both fits are printed, with the sensitivity of lambda_g to the inner limit.
Usage: python3 milky_way.py   (writes milky_way.json)"""
import numpy as np, json, os
here=os.path.dirname(os.path.abspath(__file__))
G=6.67408e-11; KPC=3.08567758128e19; MSUN=2e30
d=np.loadtxt(os.path.join(here,'milky_way_sofue2020.csv'),delimiter=',',comments='#')
X=d[:,0]*KPC; Y=d[:,1]*1000; m0=5e10*MSUN
H=np.arange(1e20/100,1e20*100,1e20/2)                       # the grid of galaxy_fit.py
def fit(x,y,m0):                                            # as in galaxy_fit.py, with the profile of the objective
    m=np.arange(m0/15,m0*15,10**np.floor(np.log10(m0))/2); best=(np.inf,None); prof=[]
    for h in H:
        mm=m[:,None]; xx=x[None,:]; ft=h/(mm/xx+h)*np.sqrt(G*mm/xx)
        w=(y[-1]-ft[:,-1])/x[-1]; dd=y[None,:]-ft-w[:,None]*xx; sd=dd.std(axis=1); k=int(np.argmin(sd))
        prof.append((h,sd[k]))
        if sd[k]<best[0]: best=(sd[k],(h,m[k],w[k],np.sqrt((dd[k]**2).mean())))
    near=[h for h,s in prof if s<=1.1*best[0]]              # lambda_g within 10 % of the minimum standard deviation
    sd,(h,mk,w,rms)=best
    return dict(points=len(x),r_min_kpc=float(x.min()/KPC),r_max_kpc=float(x.max()/KPC),lambda_g=float(h),
                lambda_g_10pc=[float(min(near)),float(max(near))],m_Msun=float(mk/MSUN),
                omega_km_s_kpc=float(w*KPC/1000),std_km_s=float(sd/1000),rms_km_s=float(rms/1000))
out={'outside the bulge, R >= 2.5 kpc':fit(X[X>=2.5*KPC],Y[X>=2.5*KPC],m0),'whole curve':fit(X,Y,m0)}
out['inner limit 2 to 5 kpc']={'%.1f kpc'%r:fit(X[X>=r*KPC*0.9999],Y[X>=r*KPC*0.9999],m0)['lambda_g'] for r in (2,2.5,3,3.5,4,5)}
print('Milky Way, Sofue (2020); mass scale 5e10 Msun; Eq. (91) by the procedure of galaxy_fit.py')
for k in ('outside the bulge, R >= 2.5 kpc','whole curve'):
    f=out[k]
    print('%-32s %2d points, %.2f to %.1f kpc: lambda_g = %.2e kg/m (%.2e to %.2e within 10 %% of the minimum),'
          ' m = %.2e Msun, omega = %.2f km/s/kpc, std %.1f km/s, rms %.1f km/s'
          %(k,f['points'],f['r_min_kpc'],f['r_max_kpc'],f['lambda_g'],*f['lambda_g_10pc'],f['m_Msun'],f['omega_km_s_kpc'],
            f['std_km_s'],f['rms_km_s']))
print('lambda_g for inner limits of 2 to 5 kpc: '+', '.join('%s %.2e'%(k,v) for k,v in out['inner limit 2 to 5 kpc'].items()))
json.dump(out,open(os.path.join(here,'milky_way.json'),'w'),indent=1)
