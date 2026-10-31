"""Spinning-disc readings (Section 5.4, Table 6): pointer positions read from the videos, converted with the
coin calibration of IMG_8699 (2.50 g = 0.2087 window widths). + = pan lighter.
Usage: python3 readings.py"""
import numpy as np
gN=9.80665
kap=2.50/(0.5693-0.3606)          # g per window width of the pointer scale (IMG_8699 coin)
r=kap/(2.50/0.329)                # transfer scale for IMG_8693 and IMG_8694 (+-25 %)
c2=2.06*1.374/154*kap
rows=[]
# runs with pointer positions in pixels of a 154 px window (IMG_8695) or in transferred units (IMG_8693, IMG_8694)
for rn,n,dm,e in [('IMG_8695',2369,2.11/154*kap,r*0.1201),('IMG_8695',2470,28.99/154*kap,r*0.2204),('IMG_8695',2591,33.29/154*kap,r*0.1797),
                  ('IMG_8694',2444,r*1.05,r*0.4034),('IMG_8694',2692,r*2.03,r*0.5173),('IMG_8693',2503,r*0.95,r*0.2510),('IMG_8693',2751,r*2.21,r*0.5596)]:
    rows.append((rn,n,dm,e))
def sube(rn,v): return np.sqrt(0.15**2+(0.03*v)**2) if rn=='IMG_8695' else max(0.35,0.25*abs(v))
for rn,n,v in [('IMG_8695',1785,15.09/154*kap),('IMG_8695',2140,-15.81/154*kap),('IMG_8694',1413,6.1*c2),('IMG_8694',2180,-4.3*c2),('IMG_8693',2062,-3.2*c2)]:
    rows.append((rn,n,v,sube(rn,v)))
rows.append(('IMG_8699',1904,0.1016*kap,0.10)); rows.append(('IMG_8699',2372,0.0328*kap,np.sqrt(0.12**2+(0.03*0.393)**2)))
print('calibration %.4f g per window width; transfer scale %.4f'%(kap,r))
print('run        n (rpm)   dm (g)            F (mN)')
for rn,n,v,e in sorted(rows,key=lambda q:q[1]): print('%s  %5d   %+.2f +- %.2f     %+.1f +- %.1f'%(rn,n,v,e,v*gN,e*gN))
print('IMG_8699   3056   >= %.2f (beam stop)  >= %.1f'%(0.194*kap,0.194*kap*gN))
print('beam stop 2.49 g = %.1f mN; 2591 rpm at 0.667 W: %.1f N/kW'%(2.49*gN,2.589e-3*gN/0.667*1e3))
