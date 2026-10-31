"""Phone magnetometer readings along the container (Section 5.3): mean field at each 10 cm hold.
Usage: python3 holds.py   (phone flat; Bx, By, Bz in the phone's axes, microtesla)"""
import numpy as np, os
here=os.path.dirname(os.path.abspath(__file__))
def load(f):
    d=np.genfromtxt(os.path.join(here,f),delimiter=',',skip_header=1,names=True); return d['time'],np.vstack([d['Bx'],d['By'],d['Bz']]).T
runs={'magnetometer_2026-10-29_00-26-09_steps_box.csv':{'0 to 60 cm on the box':[(0.2,2.5),(3.0,3.4),(5.1,7.5),(8.5,11.9),(12.7,16.5),(17.2,19.2),(19.7,21.8)]},
      'magnetometer_2026-10-29_00-38-56_steps_plank.csv':{'0 to 50 cm on the plank, reading A':[(0.75,3.7),(4.2,8.75),(9.0,13.5),(13.9,15.5),(16.0,17.7),(18.6,22.5)],
                                                            '0 to 50 cm on the plank, reading B':[(0.75,3.7),(4.2,8.75),(9.0,13.5),(13.9,17.9),(18.6,22.5),(22.9,23.22)]}}
for f,sets in runs.items():
    t,B=load(f)
    for name,segs in sets.items():
        print(f'\n{f}: {name}\n pos(cm)  |Bh|    |Bz|   |B|    turn vs 0 cm (deg)')
        B0=None
        for k,(a,b) in enumerate(segs):
            m=B[(t>=a)&(t<=b)].mean(0); B0=m if B0 is None else B0
            turn=np.degrees(np.arccos(np.clip(m@B0/np.linalg.norm(m)/np.linalg.norm(B0),-1,1)))
            print('  %3d   %6.2f  %6.2f  %6.2f   %4.1f'%(10*k,np.hypot(m[0],m[1]),-m[2],np.linalg.norm(m),turn))
t,B=load('magnetometer_2026-10-28_23-54-16_slide.csv'); n=np.linalg.norm(B,axis=1)
print('\ncontinuous slide: |B| %.2f at the start, %.2f at the end'%(n[t<1.5].mean(),n[t>11.6].mean()))
