"""Effects of the apparatus on the period difference (Section 5.2, Table 4).
Usage: python3 corrections.py   (after periods.py)"""
import json, os, numpy as np
here=os.path.dirname(os.path.abspath(__file__)); P=json.load(open(os.path.join(here,'periods.json')))
g=9.806; L=1.96; m_rope=3.0; M_steel=66.0; M_mag=184.0        # g; L = pivot to centre of the steel sphere (m)
Ts,Tm=P['IMG_8754_steel']['mean_period'],P['IMG_8755_magnets']['mean_period']
d0=100*(Tm/Ts-1); e0=100*np.hypot(P['IMG_8754_steel']['mean_error']/Ts,P['IMG_8755_magnets']['mean_error']/Tm)
rope=lambda M: np.sqrt((M+m_rope/3)/(M+m_rope/2))                 # uniform rope of mass m_rope under a bob M
c_rope=100*(rope(M_mag)/rope(M_steel)-1)           # the rope lengthens the magnets' period relative to the sphere's
c_decay,e_decay=-0.10,0.10                                       # slower decay of the magnets' swing (estimated)
r_ball=(3*(M_steel/7.85)/(4*np.pi))**(1/3)/100                   # m, 66 g steel sphere
dz=r_ball-0.01904/2                                              # flat stack 19.04 mm tall: its centre is dz higher
c_centre=+100*dz/(2*L); e_centre=0.05
A=m_rope/(1.14*200)*1e-4                                         # m^2, solid nylon cross-section
dF=(M_mag-M_steel)/1000*g*(1+np.radians(20.9)**2/4)             # extra mean tension under the magnets (N)
stretch=[100*dF/(E*A)*L/(2*L) for E in (2e9,0.5e9)]              # nylon twine, modulus 0.5 to 2 GPa
rem=d0-c_rope+c_decay+c_centre; erem=np.sqrt(e0**2+e_decay**2+e_centre**2)
print('measured                 %+.2f +- %.2f %%'%(d0,e0))
print('rope mass                %+.2f %%'%(-c_rope))
print('slower decay (estimate)  %+.2f +- %.2f %%'%(c_decay,e_decay))
print('centre %.1f mm higher     %+.2f %%'%(1000*dz,c_centre))
print('remaining                %+.2f +- %.2f %%'%(rem,erem))
print('twine stretch            %.2f to %.2f %%'%(-stretch[0],-stretch[1]))
lo,hi=rem-stretch[1],rem-stretch[0]
print('remaining with stretch   %+.2f to %+.2f %% (+- %.2f)'%(lo,hi,erem))
print('m_g/m_i of the magnets   %.4f to %.4f (+- %.4f)'%(1/(1+hi/100)**2,1/(1+lo/100)**2,2*erem/100))
