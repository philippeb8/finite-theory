"""Descriptive fits of the spinning-disc readings above the change of sign (Section 5.4, Eq. 100, Fig. 13).
dm = A [ (n/n0)^(5/2) - 1 ], weighted by the reading uncertainties; parameter errors scaled by sqrt(chi2/dof) when
chi2/dof > 1. The scaling does not correct clipping, an unknown zero or a shared calibration bias, and the rpm
uncertainty (read from the acoustic tone) is not included. The readings were selected after inspecting the data:
those of runs IMG_8693 to IMG_8695 at 2369 rpm and above, above the change of sign.
  principal fit : the four readings clear of the beam stop (2369 to 2503 rpm)
  sensitivity   : adding IMG_8695 at 2591 rpm (2.59 g, within its uncertainty of the 2.49 g stop; whether the pointer
                  was free is not established), and adding the provisional IMG_8693 / IMG_8694 readings at 2751 and
                  2692 rpm (transferred scale +-25 %, rest position and stop not known)
Usage: python3 fit_law.py   (writes fit_law.json, read by figures/make_figures.py)"""
import numpy as np, json, os
from scipy.optimize import least_squares
here=os.path.dirname(os.path.abspath(__file__))
kap=2.50/(0.5693-0.3606); r=kap/(2.50/0.329)          # coin calibration of IMG_8699 and the transfer scale, as in readings.py
pts=[('IMG_8695',2369,2.11/154*kap,r*0.1201),('IMG_8694',2444,r*1.05,r*0.4034),('IMG_8695',2470,28.99/154*kap,r*0.2204),
     ('IMG_8693',2503,r*0.95,r*0.2510),('IMG_8695',2591,33.29/154*kap,r*0.1797),('IMG_8694',2692,r*2.03,r*0.5173),
     ('IMG_8693',2751,r*2.21,r*0.5596)]
def fit(sel):
    n=np.array([p[1] for p in sel],float); y=np.array([p[2] for p in sel]); s=np.array([p[3] for p in sel]); dof=len(n)-2
    rr=least_squares(lambda p:(p[0]*((n/p[1])**2.5-1)-y)/s,[8,2330]); cov=np.linalg.inv(rr.jac.T@rr.jac); chi2=float(np.sum(rr.fun**2))
    e=np.sqrt(np.diag(cov))*max(np.sqrt(chi2/dof),1.0)
    W=1/s**2; M=np.vstack([np.ones_like(n),n]).T; c=np.linalg.solve(M.T@(M*W[:,None]),M.T@(W*y)); chil=float(np.sum((M@c-y)**2*W))
    return dict(N=len(n),A=rr.x[0],eA=e[0],n0=rr.x[1],en0=e[1],chi2=chi2,dof=dof,line_zero=-c[0]/c[1],line_chi2=chil)
out={'principal, 4 readings clear of the stop':fit([p for p in pts if p[1]<=2503]),
     'with 2591 rpm, 5 readings':fit([p for p in pts if p[1]<=2591]),
     'with 2591, 2692 and 2751 rpm, 7 readings':fit(pts)}
for k,v in out.items():
    print('%-42s A = %5.2f +- %4.2f g   n0 = %4.0f +- %2.0f rpm   chi2 = %5.2f / %d   line through zero at %4.0f rpm (chi2 %.2f)'
          %(k,v['A'],v['eA'],v['n0'],v['en0'],v['chi2'],v['dof'],v['line_zero'],v['line_chi2']))
p=out['principal, 4 readings clear of the stop']
print('principal fit at 2372 rpm: %.2f g (IMG_8699 measured 0.39 +- 0.12 g, not used in the fit)'%(p['A']*((2372/p['n0'])**2.5-1)))
json.dump(out,open(os.path.join(here,'fit_law.json'),'w'),indent=1)
