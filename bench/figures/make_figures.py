"""Figures 9, 10, 12 and 13 of the paper from the records in this folder.
Usage: python3 make_figures.py   (after freefall/ff_refit.py, pendulum/periods.py and pendulum/corrections.py; writes out/*.eps|pdf|png)"""
import numpy as np, json, glob, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import os, glob as _g
for f in _g.glob('/usr/share/fonts/truetype/crosextra/Carlito-*.ttf'): fm.fontManager.addfont(f)
HAVE_CARLITO=bool(_g.glob('/usr/share/fonts/truetype/crosextra/Carlito-*.ttf'))
plt.rcParams.update({'font.family':'Carlito' if HAVE_CARLITO else 'DejaVu Sans','font.size':9,'axes.edgecolor':'#9aa3b2','axes.labelcolor':'#0C1A33',
 'xtick.color':'#4a5568','ytick.color':'#4a5568',
 'ps.fonttype':42,'pdf.fonttype':42,'axes.linewidth':0.8})
INK='#0C1A33'; CYAN='#137F99'; AMBER='#C56A00'; MUTE='#6b7280'; GRID='#e6e9ef'; LIGHT='#b8bfcc'
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.join(HERE,'..')
OUT=os.path.join(HERE,'out')+'/'; os.makedirs(OUT,exist_ok=True)
def style(a):
    a.grid(color=GRID,lw=0.7); a.set_axisbelow(True)
    for s in ('top','right'): a.spines[s].set_visible(False)
def save(fig,name):
    for ext in ('eps','pdf','png'): fig.savefig(OUT+name+'.'+ext,dpi=200 if ext=='png' else None,bbox_inches='tight')

# ---------------- 1. Free fall, no prediction bands ----------------
J=json.load(open(os.path.join(ROOT,'freefall','lower_half.json')))
B=J['blocks']; geff=J['geff']
fig,(a,b)=plt.subplots(1,2,figsize=(6.6,3.2),gridspec_kw={'width_ratios':[1.55,1]})
labs=[f'{d}{o}' for d in (7,8,9,10,11) for o in 'UD']
for k,lab in enumerate(labs):
    U=lab.endswith('U'); col=AMBER if U else CYAN; lt='#E9C49C' if U else '#A3CDD7'; mk='o' if U else 's'
    v=np.array(B[lab]); a.plot(k+np.linspace(-0.15,0.15,len(v)),v,mk,ms=3.6,color=lt,mec='none',zorder=2)
    a.errorbar([k+0.32],[v.mean()],yerr=[v.std(ddof=1)/np.sqrt(len(v))],fmt=mk,ms=5.5,color=col,mec='white',mew=0.8,ecolor=col,elinewidth=1.2,capsize=2,zorder=3)
a.axhline(geff,color=MUTE,lw=1,ls='--',zorder=1); a.text(-0.5,geff-0.3,'$g$ less air drag',fontsize=7.4,color=MUTE,ha='left',va='top')
a.plot([],[],'o',color=AMBER,ms=5,label='U drops'); a.plot([],[],'s',color=CYAN,ms=5,label='D drops')
a.legend(loc='lower left',fontsize=7.6,frameon=False,ncol=2)
a.set_xticks(range(10)); a.set_xticklabels(labs); a.set_xlabel('Session block (day, orientation U/D)')
a.set_ylabel('Magnet acceleration, lower half (m/s$^2$)'); a.set_ylim(0,14.5); a.set_xlim(-0.6,9.7); style(a)
a.set_title('(a) Lower half of each drop, calibrated',loc='left',fontsize=9,color=INK)
odd=np.array(J['odd']); even=np.array(J['even'])
b.axhline(0,color=MUTE,lw=0.9,zorder=1)
jit=np.linspace(-0.12,0.12,5)
b.plot(0+jit,even,'o',ms=5,color='#8a93a3',mec='white',mew=0.6,zorder=3)
b.plot(1+jit,odd,'o',ms=5,color='#8a93a3',mec='white',mew=0.6,zorder=3)
for i,v in enumerate((even,odd)):
    m=v.mean(); se=v.std(ddof=1)/np.sqrt(len(v))
    b.errorbar([i+0.3],[m],yerr=[se],fmt='D',ms=5.8,color=AMBER,mec='white',mew=0.8,ecolor=AMBER,elinewidth=1.4,capsize=3,zorder=4)
    b.text(i+0.42,m,f'{m:+.2f} ± {se:.2f}'.replace('-','−'),va='center',fontsize=8,color=INK)
b.set_xticks([0,1]); b.set_xticklabels(['even\n(U+D)/2 − $g$','odd\n(D − U)/2']); b.set_xlim(-0.45,1.95)
b.set_ylabel('Component (m/s$^2$)'); b.set_ylim(-4.5,1.6); style(b)
b.set_title('(b) Orientation components',loc='left',fontsize=9,color=INK)
plt.tight_layout(w_pad=1.5); save(fig,'FreeFall_Components'); plt.close(fig)
print('ff even %.2f odd %.2f'%(even.mean(),odd.mean()))

# ---------------- 2. Spinning disc, readings only ----------------
gN=9.80665
kap=2.50/(0.5693-0.3606); r=kap/(2.50/0.329); c2=2.06*1.374/154*kap
n=np.array([2369,2470,2591,2444,2692,2503,2751.]); y=np.r_[np.array([2.11,28.99,33.29])/154*kap, r*np.array([1.05,2.03,0.95,2.21])]
s=r*np.array([0.1201,0.2204,0.1797,0.4034,0.5173,0.2510,0.5596]); sn=np.array([8,8,8,10,10,10,10.])
runs=['IMG_8695']*3+['IMG_8694']*2+['IMG_8693']*2
sub=[('IMG_8695',1785,15.09/154*kap),('IMG_8695',2140,-15.81/154*kap),('IMG_8694',1413,6.1*c2),('IMG_8694',2180,-4.3*c2),('IMG_8693',2062,-3.2*c2)]
def suberr(rn,v): return np.sqrt(0.15**2+(0.03*v)**2) if rn=='IMG_8695' else max(0.35,0.25*abs(v))
mk={'IMG_8695':'o','IMG_8694':'s','IMG_8693':'^'}; col={'IMG_8695':INK,'IMG_8694':CYAN,'IMG_8693':'#4a5568'}
fig,a=plt.subplots(figsize=(6.2,3.6))
a.axhline(0,color=MUTE,lw=0.8,zorder=1)
a.axhline(2.49*gN,color=MUTE,lw=1,ls='--',zorder=1); a.text(1365,2.49*gN-2.6,'beam stop: largest measurable lightening, 24.4 mN',fontsize=7.6,color=MUTE)
FL=json.load(open(os.path.join(ROOT,'spinning_disc','fit_law.json')))
_p=FL['principal, 4 readings clear of the stop']; _q=FL['with 2591, 2692 and 2751 rpm, 7 readings']
_nn=np.linspace(_p['n0'],2520,120)
a.plot(_nn,_p['A']*((_nn/_p['n0'])**2.5-1)*gN,'--',color=CYAN,lw=1.4,zorder=2,label='descriptive fit, 4 readings clear of the stop')
_nn=np.linspace(_q['n0'],2780,120)
a.plot(_nn,_q['A']*((_nn/_q['n0'])**2.5-1)*gN,':',color=CYAN,lw=1.1,zorder=2,label='sensitivity fit with the 3 uncertain readings')
for rn in ('IMG_8695','IMG_8694','IMG_8693'):
    ii=[i for i,q in enumerate(runs) if q==rn]
    nv=np.r_[n[ii],[q[1] for q in sub if q[0]==rn]]; yv=np.r_[y[ii],[q[2] for q in sub if q[0]==rn]]
    ev=np.r_[s[ii],[suberr(rn,q[2]) for q in sub if q[0]==rn]]; o=np.argsort(nv)
    a.errorbar(nv[o],yv[o]*gN,yerr=ev[o]*gN,fmt=mk[rn]+'-',lw=0.7,ms=5.2,color=col[rn],mec='white',mew=0.8,ecolor=LIGHT,elinewidth=1,capsize=2,zorder=4,label=rn)
    for nq,yq in zip(nv[o],yv[o]):                                     # readings whose point value is not established
        if nq in (2591,2692,2751): a.plot([nq],[yq*gN],mk[rn],ms=5.2,mfc='white',mec=col[rn],mew=1.1,zorder=5)
a.errorbar([1904,2372],[0.1016*kap*gN,0.0328*kap*gN],yerr=[0.10*gN,np.sqrt(0.12**2+(0.03*0.393)**2)*gN],fmt='D-',lw=0.7,ms=5.8,color=AMBER,mec='white',mew=0.8,ecolor=AMBER,capsize=2.5,zorder=5,label='IMG_8699 (calibrated in frame)')
a.plot([2372,3056],[0.0328*kap*gN,0.194*kap*gN],':',color=AMBER,lw=0.8)
a.plot([3056],[0.194*kap*gN],'D',ms=5.8,mfc='white',mec=AMBER,mew=1.2,zorder=5)
a.annotate('',xy=(3056,0.194*kap*gN+11),xytext=(3056,0.194*kap*gN+0.9),arrowprops=dict(arrowstyle='->',color=AMBER,lw=1.4))
a.text(3150,0.194*kap*gN-9.5,'≥ 22.8 mN\n(held at stop)',ha='right',fontsize=8,color=AMBER)
a.set_xlim(1350,3180); a.set_ylim(-16,52)
a.set_xlabel('Disc speed $n$ (rpm)'); a.set_ylabel('Apparent axial force (mN), + = pan lighter')
a.plot([],[],'o',ms=5.2,mfc='white',mec=MUTE,mew=1.1,label='open: near or beyond the stop, point value not established')
a.legend(loc='upper left',ncol=2,fontsize=6.6,frameon=False,columnspacing=1.0,handlelength=2.2); style(a)
plt.tight_layout(); save(fig,'CD_Axial_Force'); plt.close(fig)

# ---------------- 3. Pendulum ----------------
P=json.load(open(os.path.join(ROOT,'pendulum','periods.json')))
fig,(a,b)=plt.subplots(1,2,figsize=(6.6,3.0),gridspec_kw={'width_ratios':[1.45,1]})
C=json.load(open(os.path.join(ROOT,'pendulum','corrections.json')))
Tm={n_:P[n_]['mean_period'] for n_ in ('IMG_8754_steel','IMG_8755_magnets')}
for side,mk,lab_ in (('near','o','near end of the swing'),('far','s','far end of the swing')):
    ts=np.array(P['IMG_8754_steel'][side]['turning_times']); tmg=np.array(P['IMG_8755_magnets'][side]['turning_times'])
    k=np.arange(len(ts)); lag=(tmg-tmg[0])-(ts-ts[0])           # extra time the magnets take, from the first turning point
    a.plot(k,lag,mk,ms=4.0 if mk=='o' else 3.6,color=AMBER if side=='near' else CYAN,mec='white',mew=0.5,zorder=3,label='Turning points, '+lab_)
dT=Tm['IMG_8755_magnets']-Tm['IMG_8754_steel']
a.plot([0,20],[0,20*dT],color=INK,lw=1.1,zorder=2,label='Difference of the periods, %.4f s per swing'%dT)
a.set_xlim(-0.6,21.2); a.set_ylim(-0.12,0.66); a.set_xticks([0,5,10,15,20])
a.set_xlabel('Swings after the first turning point'); a.set_ylabel('Extra time taken by the magnets (s)'); style(a)
a.legend(loc='upper left',fontsize=7.2,frameon=False)
a.text(20.8,0.02,'same 41 turning points\nin both records',fontsize=7.4,color=MUTE,ha='right',va='bottom')
a.set_title('(a) Same rope, %.2f m, released at about 21°'%C['length_m'],loc='left',fontsize=9,color=INK)
lab=['measured','rope mass\nremoved','swing decay\nremoved','centre height\nremoved','twine stretch\nremoved']
val=[C[k][0] for k in ('measured','after_rope','after_decay','after_centre')]
err=[C[k][1] for k in ('measured','after_rope','after_decay','after_centre')]
lo,hi,ew=C['with_stretch']
b.barh([4,3,2,1],val,color=INK,height=0.62,xerr=err,ecolor=MUTE,capsize=2,error_kw=dict(lw=0.9))
b.barh([0],[hi-lo],left=[lo],color=AMBER,height=0.62,xerr=[[0],[ew]],ecolor=MUTE,capsize=2,error_kw=dict(lw=0.9))
b.set_xlim(0,1.0)
for yy,v,e_ in zip([4,3,2,1],val,err): b.text(v+e_+0.03,yy,f'{v:.2f} %',va='center',fontsize=7.8,color=INK)
b.text(hi+ew+0.03,0,'%.2f to %.2f %%'%(lo,hi),va='center',fontsize=7.8,color=AMBER)
b.set_yticks([4,3,2,1,0]); b.set_yticklabels(lab,fontsize=7.6)
b.set_xlabel("Magnets' period longer than steel's (%)"); style(b)
b.set_title('(b) Period difference',loc='left',fontsize=9,color=INK)
plt.tight_layout(w_pad=1.2); save(fig,'Pendulum_Periods'); plt.close(fig)
print('pendulum periods', Tm)

# ---------------- 4. Float at night and on the plastic support ----------------
fig,(a,b)=plt.subplots(1,2,figsize=(6.6,3.0),gridspec_kw={'width_ratios':[1.2,1]})
PX=248/15.0
from scipy.ndimage import median_filter
sty={'IMG_8745_Double-South':(AMBER,'o','8745, two blocks, S'),'IMG_8746_Double-North':(AMBER,'s','8746, two blocks, N'),
     'IMG_8747_Single-South':(CYAN,'o','8747, one block, S'),'IMG_8748_Single-North':(CYAN,'s','8748, one block, N')}
for f in sorted(glob.glob(os.path.join(ROOT,'floating_magnet','tracks','IMG_874*_track.csv'))):
    nm=f.split('/')[-1].replace('_track.csv',''); A=np.loadtxt(f,delimiter=','); t,x,w=A[:,0],A[:,1],A[:,2]
    ok=np.abs(w-np.nanmedian(w))<15; ok&=np.abs(x-np.interp(t,t[ok],x[ok]))<8
    t,x=t[ok],x[ok]; d=(x-x[:10].mean())/PX
    c_,m_,l_=sty[nm]; sel=np.arange(0,len(t),12)
    a.plot(t[sel],d[sel],m_,ms=3.0,color=c_,mec='white',mew=0.4,label=l_)
a.set_xlabel('Time (s)'); a.set_ylabel('Displacement toward south (cm)'); style(a)
a.legend(loc='upper left',fontsize=7.2,frameon=False)
a.set_title('(a) 22 h, container on plastic boxes',loc='left',fontsize=9,color=INK)
PXP=990/60.0                                                  # plastic support: 60 cm container over 990 pixels
P1=np.loadtxt(os.path.join(ROOT,'floating_magnet','tracks','IMG_8765_plastic_track.csv'),delimiter=',')
P2=np.loadtxt(os.path.join(ROOT,'floating_magnet','tracks','IMG_8766_plastic_track.csv'),delimiter=',')
P2=P2[P2[:,0]<=98]                                            # to the end of the fitted range (fit_drift.py), before the wall
b.axvspan(0,15,color=GRID,lw=0,zorder=0); b.text(7.5,-3.2,'turn of\nthe two\nblocks',fontsize=7.2,color=MUTE,ha='center',va='top')
b.axhline(0,color=MUTE,lw=0.8,ls='--')
b.plot(P1[:,0],(P1[:,1]-P1[0,1])/PXP,'-',color=CYAN,lw=1.5,label='8765, one block: west')
b.plot(P2[:,0],(P2[:,1]-P2[0,1])/PXP,'-',color=AMBER,lw=1.5,label='8766, two blocks: east')    # left edge of the float
b.set_xlabel('Time (s)'); b.set_ylabel('Displacement toward west (cm)'); style(b)
b.legend(loc='upper left',fontsize=7.2,frameon=False)
b.set_title('(b) Plastic support, no steel below',loc='left',fontsize=9,color=INK)
plt.tight_layout(w_pad=1.2); save(fig,'Floating_Magnet_Night'); plt.close(fig)
print('float ok')
