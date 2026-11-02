"""Refit of the 40 paired free-fall drops (Section 5.1, Table 3 and Fig. 9 of the paper).
Usage: python3 ff_refit.py data   (folder holding day7..day11/FreeFall-7pm-{Up,Down}.py)
Reads the DATA blocks of the ten scripts (does not execute them)."""
import sys, os, re, numpy as np
from scipy import stats
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
SRC=sys.argv[1] if len(sys.argv)>1 else '.'; TS=30/240; g=9.806; geff=g*(1-0.0036)
def parse(lines):
    tr=[[] for _ in range(4)]
    for l in lines:
        tok=l.replace('\t',' ').split()
        for i in range(4):
            try: tr[i].append((float(tok[3*i])*TS,float(tok[3*i+1]),float(tok[3*i+2])))
            except (ValueError,IndexError): pass
    return [np.array(t) for t in tr if len(t)>=5]
drops=[]
for d in (7,8,9,10,11):
    for o,fn in (('U','FreeFall-7pm-Up.py'),('D','FreeFall-7pm-Down.py')):
        L=[l for l in re.search(r'DATA = r"""(.*?)"""',open(os.path.join(SRC,f'day{d}',fn)).read(),re.S).group(1).split('\n') if l.strip()]
        iS=L.index('Steel')
        for m,s in zip(parse(L[2:iS]),parse(L[iS+2:])): drops.append(dict(d=d,o=o,m=m,s=s))
def quad(t,y):
    A=np.vstack([np.ones_like(t),t,-0.5*t**2]).T; c=np.linalg.lstsq(A,y,rcond=None)[0]; return c,y-A@c
for r in drops:
    t=r['m'][:,0]; r['t']=t; r['aM']=quad(t,r['m'][:,2])[0][2]; r['aS']=quad(t,r['s'][:,2])[0][2]; r['R']=r['aM']/r['aS']
bad=[r for r in drops if r['aS']<9]; drops=[r for r in drops if r['aS']>=9]
print('excluded (failed steel track):',[(r['d'],r['o'],round(r['aS'],2)) for r in bad])
# ---- ratio method (paper) ----
print('\nblock  n   a_magnet        a_steel         ratio')
blk={}
for d in (7,8,9,10,11):
    for o in 'UD':
        b=[r for r in drops if r['d']==d and r['o']==o]; n=len(b); blk[(d,o)]=np.mean([r['R'] for r in b])
        f=lambda k: (np.mean([r[k] for r in b]),np.std([r[k] for r in b],ddof=1)/np.sqrt(n))
        print('%2d%s   %d   %.2f ± %.2f   %.2f ± %.2f   %.3f ± %.3f'%(d,o,n,*f('aM'),*f('aS'),*f('R')))
Rr=np.array([r['R'] for r in drops]); N=len(Rr); S=np.array([r['aS'] for r in drops])
print('\nsteel %.3f ± %.3f m/s2 (%.1f%% above g)'%(S.mean(),S.std(ddof=1)/np.sqrt(N),100*(S.mean()/g-1)))
print('ratio %.4f ± %.4f ; deficit %.3f ± %.3f m/s2 (%.1f ± %.1f%%)'%(Rr.mean(),Rr.std(ddof=1)/np.sqrt(N),g*(1-Rr.mean()),g*Rr.std(ddof=1)/np.sqrt(N),100*(1-Rr.mean()),100*Rr.std(ddof=1)/np.sqrt(N)))
U=np.array([r['R'] for r in drops if r['o']=='U']); D=np.array([r['R'] for r in drops if r['o']=='D'])
print('sd of ratio: U %.3f, D %.3f ; Levene p = %.4f'%(U.std(ddof=1),D.std(ddof=1),stats.levene(U,D).pvalue))
odd=np.array([g*(blk[(d,'D')]-blk[(d,'U')])/2 for d in (7,8,9,10,11)]); even=np.array([g*((blk[(d,'D')]+blk[(d,'U')])/2-1) for d in (7,8,9,10,11)])
tt=stats.ttest_1samp(odd,0)
print('odd per day',np.round(odd,3),' mean %.3f ± %.3f m/s2, t = %.2f, p = %.3f'%(odd.mean(),odd.std(ddof=1)/np.sqrt(5),tt.statistic,tt.pvalue))
print('even mean %.3f ± %.3f m/s2'%(even.mean(),even.std(ddof=1)/np.sqrt(5)))
# ---- optics calibrated on the steel sphere: y_obs - y_ref = s (Z + k Z^2) ----
Zf=lambda t,Z0,V0: Z0+V0*t-0.5*geff*t**2
for r in drops:
    t=r['t']; yS=r['s'][:,2]; yM=r['m'][:,2]; yref=yS[0]
    sol=least_squares(lambda p: p[0]*(Zf(t,p[2],p[3])+p[1]*Zf(t,p[2],p[3])**2)-(yS-yref),[1.08,0.0,0.0,-0.3])
    s_,k_=sol.x[:2]; w=(yM-yref)/s_; r['Zm']=(-1+np.sqrt(1+4*k_*w))/(2*k_); r['k']=k_
def law(r,name):
    t,Z=r['t'],r['Zm']
    if name=='const': c,res=quad(t,Z); return res@res,c[2]
    if name=='lin':
        f=lambda p: p[0]-geff/p[2]*t+(p[1]+geff/p[2])/p[2]*(1-np.exp(-p[2]*t))-Z
        s=least_squares(f,[Z[0],-0.3,0.4],bounds=([-1,-3,1e-4],[1,3,20])); return s.fun@s.fun,s.x[2]
    rhs={'quad':lambda p:(lambda tt,u:[u[1],-geff+p[2]*u[1]**2]),'tn4':lambda p:(lambda tt,u:[u[1],-geff+p[2]*tt*u[1]**2])}[name]
    f=lambda p: p[0]+solve_ivp(rhs(p),[0,t[-1]],[0,p[1]],t_eval=t,rtol=1e-9,atol=1e-11).y[0]-Z
    s=least_squares(f,[Z[0],-0.3,0.1],bounds=([-1,-3,0],[1,3,5])); return s.fun@s.fun,s.x[2]
print('\noptics curvature k = %.4f ± %.4f /m'%(np.mean([r['k'] for r in drops]),np.std([r['k'] for r in drops])))
res={}
for name in ('lin','quad','const','tn4'):
    out=[law(r,name) for r in drops]; res[name]=(sum(o[0] for o in out),out)
npts=sum(len(r['t']) for r in drops); s2=res['lin'][0]/(npts-3*N)
for name,(ssr,out) in res.items():
    pu=[o[1] for o,r in zip(out,drops) if r['o']=='U']; pd=[o[1] for o,r in zip(out,drops) if r['o']=='D']
    print('%-5s dchi2 %+5.0f   U %.3f ± %.3f   D %.3f ± %.3f'%(name,(ssr-res['lin'][0])/s2,np.mean(pu),np.std(pu,ddof=1)/np.sqrt(len(pu)),np.mean(pd),np.std(pd,ddof=1)/np.sqrt(len(pd))))
# ---- lateral motion of the magnet relative to the steel sphere's drift at the same height ----
print()
for o in 'UD':
    acc=[];end=[]
    for r in [r for r in drops if r['o']==o]:
        t=r['t']; xM,yM=r['m'][:,1],r['m'][:,2]; xS,yS=r['s'][:,1],r['s'][:,2]; i=np.argsort(yS)
        dx=(xM-xM[0]*np.interp(yM,yS[i],(xS/xS[0])[i]))*np.sign(-xM[0])   # + toward the steel sphere
        A=np.vstack([np.ones_like(t),t,0.5*t**2]).T; acc.append(np.linalg.lstsq(A,dx,rcond=None)[0][2]); end.append(dx[-1])
    print('lateral %s: acceleration %.2f ± %.2f m/s2, end offset %.3f ± %.3f m'%(o,np.mean(acc),np.std(acc,ddof=1)/np.sqrt(len(acc)),np.mean(end),np.std(end,ddof=1)/np.sqrt(len(end))))

# ---- lower half only (a wall stands beside the upper half): two-half fit, continuous at mid-height ----
def knot(t,Z,tc,sign=-1):
    A=np.vstack([np.ones_like(t),t,sign*0.5*t**2,sign*0.5*np.clip(t-tc,0,None)**2]).T
    c=np.linalg.lstsq(A,Z,rcond=None)[0]; return c[2],c[2]+c[3]
def cross(t,Z,zc):
    i=np.argmax(Z<=zc); return np.interp(zc,[Z[i],Z[i-1]],[t[i],t[i-1]])
print('\nLOWER HALF')
H={'A':[],'B':[],'R':[]}
for r in drops:
    t=r['t']; yS=r['s'][:,2]; yM=r['m'][:,2]; yref=yS[0]
    sA=least_squares(lambda p: p[0]*(Zf(t,p[2],p[3])+p[1]*Zf(t,p[2],p[3])**2)-(yS-yref),[1.08,0.0,0.0,-0.3])
    s_,k_,Z0,V0=sA.x; ZA=(-1+np.sqrt(1+4*k_*(yM-yref)/s_))/(2*k_); zc=Z0+0.5*(Zf(t[-1],Z0,V0)-Z0)
    tc=cross(t,ZA,zc); a1,a2=knot(t,ZA,tc)
    xM,xS=r['m'][:,1],r['s'][:,1]; ts=cross(t,yS,yS[0]+0.5*(yS[-1]-yS[0])); tm=cross(t,yM,yS[0]+0.5*(yS[-1]-yS[0]))
    lm=knot(t,np.sign(-xM[0])*xM,tm,+1)[1]; ls=knot(t,np.sign(-xS[0])*xS,ts,+1)[1]   # + toward the other sphere
    Ak=np.vstack([np.ones_like(t),t,-0.5*t**2,-0.5*np.clip(t-tc,0,None)**2]).T; ck=np.linalg.lstsq(Ak,ZA,rcond=None)[0]
    H['A'].append(dict(d=r['d'],o=r['o'],a1=a1,a2=a2,lat=lm-ls,s=s_,k=k_,cal_rms=float(np.sqrt(np.mean(sA.fun**2))),fit_rms=float(np.sqrt(np.mean((ZA-Ak@ck)**2))),tc=tc,n=len(t)))
    sB=least_squares(lambda p: p[0]*Zf(t+p[1]*t**3,p[2],p[3])-(yS-yref),[1.0,0.0,0.0,-0.3])
    sb,cb,Z0b,V0b=sB.x; tau=t+cb*t**3; ZB=(yM-yref)/sb; zcb=Z0b+0.5*(Zf(tau[-1],Z0b,V0b)-Z0b)
    b1,b2=knot(tau,ZB,cross(tau,ZB,zcb)); H['B'].append(dict(d=r['d'],o=r['o'],a1=b1,a2=b2))
    yc=yS[0]+0.5*(yS[-1]-yS[0]); m1,m2=knot(t,yM,cross(t,yM,yc)); s1,s2=knot(t,yS,cross(t,yS,yc))
    H['R'].append(dict(d=r['d'],o=r['o'],a1=geff*m1/s1,a2=geff*m2/s2))
for key,lab in (('A','height-calibrated (paper)'),('B','clock-calibrated'),('R','ratio, same lower half')):
    for part in ('a1','a2'):
        a=np.array([x[part] for x in H[key]]); n=len(a)
        B={(d,o):np.mean([x[part] for x in H[key] if x['d']==d and x['o']==o]) for d in (7,8,9,10,11) for o in 'UD'}
        odd=np.array([(B[(d,'D')]-B[(d,'U')])/2 for d in (7,8,9,10,11)]); even=np.array([(B[(d,'D')]+B[(d,'U')])/2-geff for d in (7,8,9,10,11)])
        print('%-26s %s: a %.2f ± %.2f (deficit %.1f ± %.1f%%) | odd %+.3f ± %.3f | even %+.2f ± %.2f'%(lab,{'a1':'upper','a2':'lower'}[part],a.mean(),a.std(ddof=1)/np.sqrt(n),100*(1-a.mean()/geff),100*a.std(ddof=1)/np.sqrt(n)/geff,odd.mean(),odd.std(ddof=1)/np.sqrt(5),even.mean(),even.std(ddof=1)/np.sqrt(5)))
for o in 'UD':
    v=np.array([x['lat'] for x in H['A'] if x['o']==o]); print('lower-half lateral acceleration of the magnet relative to the steel sphere, %s: %.2f ± %.2f m/s2 (+ toward it)'%(o,v.mean(),v.std(ddof=1)/np.sqrt(len(v))))
# ---- Table 3 of the paper (height-calibrated, per session block) and the data of Fig. 9
import json
print('\nTABLE 3  block  drops  a1 upper (m/s2)   a2 lower (m/s2)')
blocks={}
for d in (7,8,9,10,11):
    for o in 'UD':
        b=[x for x in H['A'] if x['d']==d and x['o']==o]; n=len(b)
        a1=np.array([x['a1'] for x in b]); a2=np.array([x['a2'] for x in b]); blocks[f'{d}{o}']=a2.tolist()
        print('        %2d%s     %d     %.2f +- %.2f      %.2f +- %.2f'%(d,o,n,a1.mean(),a1.std(ddof=1)/np.sqrt(n),a2.mean(),a2.std(ddof=1)/np.sqrt(n)))
B2={k:np.mean(v) for k,v in blocks.items()}
odd=[(B2[f'{d}D']-B2[f'{d}U'])/2 for d in (7,8,9,10,11)]; even=[(B2[f'{d}D']+B2[f'{d}U'])/2-geff for d in (7,8,9,10,11)]
out=os.path.join(os.path.dirname(os.path.abspath(__file__)),'lower_half.json')
json.dump({'blocks':blocks,'odd':odd,'even':even,'geff':geff},open(out,'w'),indent=1); print('wrote',os.path.relpath(out))

# ---- per-drop calibration and fit record (scale s, curvature k of the steel calibration y = s(Z + kZ^2),
#      rms of the calibration fit and of the two-half fit, mid-height time tc, number of tracked positions)
import csv
outc=os.path.join(os.path.dirname(os.path.abspath(__file__)),'per_drop.csv')
with open(outc,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['day','orientation','positions','scale_s','curvature_k_per_m','calibration_rms_m','two_half_fit_rms_m','t_mid_s','a1_upper','a2_lower','lateral_rel_steel'])
    for x in H['A']: w.writerow([x['d'],x['o'],x['n'],'%.5f'%x['s'],'%.4f'%x['k'],'%.3f'%x['cal_rms'],'%.5f'%x['fit_rms'],'%.4f'%x['tc'],'%.3f'%x['a1'],'%.3f'%x['a2'],'%.3f'%x['lat']])
print('wrote',os.path.relpath(outc))
