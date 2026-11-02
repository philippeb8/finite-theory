"""Effects of the apparatus on the period difference (Section 5.2, Table 4). The rope and centre corrections are
calculated from the geometry, the decay correction is measured, and the twine stretch is a sensitivity interval for
nylon moduli of 0.5 to 2 GPa, not a confidence interval.
Usage: python3 corrections.py   (after periods.py; writes corrections.json, read by figures/make_figures.py)"""
import json, os, numpy as np
here=os.path.dirname(os.path.abspath(__file__)); P=json.load(open(os.path.join(here,'periods.json')))
g=9.806; m_rope=3.0; M_steel=66.0; M_mag=184.0                     # g
S,M=P['IMG_8754_steel'],P['IMG_8755_magnets']
Ts,Tm=S['mean_period'],M['mean_period']
d0=100*(Tm/Ts-1); e0=100*np.hypot(S['mean_error']/Ts,M['mean_error']/Tm)
rope=lambda Mb: np.sqrt((Mb+m_rope/3)/(Mb+m_rope/2))               # uniform rope of mass m_rope under a bob Mb
L=g*(S['period_zero_amplitude']/(2*np.pi))**2/rope(M_steel)**2      # m, pivot to the centre of the steel sphere
c_rope=100*(rope(M_mag)/rope(M_steel)-1)           # the rope lengthens the magnets' period relative to the sphere's
# the swing of the magnets decays more slowly (20.5 -> 14.7 deg against 12.0 deg for the sphere): measured from the
# peak-to-peak height of the bob in the image (periods.py); +-0.02 % allows 10 % on the amplitude scale
c_decay=100*(M['period_zero_amplitude']/S['period_zero_amplitude']-Tm/Ts); e_decay=0.02
r_ball=(3*(M_steel/7.85)/(4*np.pi))**(1/3)/100                   # m, 66 g steel sphere
dz=r_ball-0.01904/2                                              # flat stack 19.04 mm tall: its centre is dz higher
c_centre=+100*dz/(2*L); e_centre=0.05                            # +-0.05 % allows +-2 mm on the length
A=m_rope/(1.14*200)*1e-4                                         # m^2, solid nylon cross-section
dF=(M_mag-M_steel)/1000*g*(1+np.radians(20.9)**2/4)             # extra mean tension under the magnets (N)
stretch=[100*dF/(E*A)/2 for E in (2e9,0.5e9)]                    # nylon twine, modulus 0.5 to 2 GPa
rem=d0-c_rope+c_decay+c_centre; erem=np.sqrt(e0**2+e_decay**2+e_centre**2)
lo,hi=rem-stretch[1],rem-stretch[0]
print('length, pivot to the centre of the sphere  %.2f m'%L)
print('measured                 %+.2f +- %.2f %%'%(d0,e0))
print('rope mass                %+.2f %%'%(-c_rope))
print('slower decay (measured)  %+.2f +- %.2f %%'%(c_decay,e_decay))
print('centre %.1f mm higher     %+.2f +- %.2f %%'%(1000*dz,c_centre,e_centre))
print('remaining                %+.2f +- %.2f %%'%(rem,erem))
print('twine stretch            %.2f to %.2f %%'%(-stretch[0],-stretch[1]))
print('remaining with stretch   %+.2f to %+.2f %% (+- %.2f)'%(lo,hi,erem))
# For an ideal point pendulum m_g/m_i = (L_mag/L_steel)(T_steel/T_mag)^2. The two readings below are conditional:
# the mass-ratio equivalent assumes equal effective lengths and unchanged non-gravitational torques; the length
# equivalent assumes an unchanged mass ratio. Neither is a measured property of the apparatus.
dL=[1000*L*((1+x/100)**2-1) for x in (lo,hi)]
print('conditional m_g/m_i, equal effective lengths   %.3f to %.3f (+- %.3f)'%(1/(1+hi/100)**2,1/(1+lo/100)**2,2*erem/100))
print('or, with m_g/m_i = 1, magnets\' effective length longer by %.1f to %.1f mm'%tuple(dL))
json.dump(dict(length_m=L,measured=[d0,e0],after_rope=[d0-c_rope,e0],after_decay=[d0-c_rope+c_decay,np.hypot(e0,e_decay)],
               after_centre=[rem,erem],with_stretch=[lo,hi,erem],mg_mi=[1/(1+hi/100)**2,1/(1+lo/100)**2,2*erem/100],
               length_equivalent_mm=dL),
          open(os.path.join(here,'corrections.json'),'w'),indent=1)
