import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# =========================================================
FPS_RECORDED = 240.0
FPS_PLAYBACK = 30.0
TIME_SCALE = FPS_PLAYBACK / FPS_RECORDED   # 0.125
# =========================================================

DATA = r"""Magnet
t	x	y	t	x	y	t	x	y	t	x	y
0.00E+00	-2.67E-01	-3.71E-02	0.00E+00	-2.67E-01	-4.45E-02	0.00E+00	-2.67E-01	-2.97E-02	0.00E+00	-2.60E-01	-4.45E-02
5.00E-01	-2.60E-01	-1.04E-01	5.00E-01	-2.67E-01	-1.19E-01	5.00E-01	-2.60E-01	-8.90E-02	5.00E-01	-2.52E-01	-1.26E-01
1.00E+00	-2.60E-01	-2.00E-01	1.00E+00	-2.67E-01	-2.30E-01	1.00E+00	-2.60E-01	-1.78E-01	1.00E+00	-2.52E-01	-2.22E-01
1.50E+00	-2.52E-01	-3.41E-01	1.50E+00	-2.60E-01	-3.63E-01	1.50E+00	-2.52E-01	-3.19E-01	1.50E+00	-2.37E-01	-3.86E-01
2.00E+00	-2.52E-01	-5.12E-01	2.00E+00	-2.60E-01	-5.34E-01	2.00E+00	-2.67E-01	-4.89E-01	2.00E+00	-2.52E-01	-5.27E-01
2.50E+00	-2.45E-01	-7.19E-01	2.50E+00	-2.60E-01	-7.42E-01	2.50E+00	-2.60E-01	-6.67E-01	2.50E+00	-2.52E-01	-7.49E-01
3.00E+00	-2.37E-01	-9.49E-01	3.00E+00	-2.60E-01	-9.79E-01	3.00E+00	-2.52E-01	-9.12E-01	3.00E+00	-2.45E-01	-9.86E-01
3.50E+00	-2.37E-01	-1.20E+00	3.50E+00	-2.52E-01	-1.26E+00	3.50E+00	-2.45E-01	-1.18E+00	3.50E+00	-2.45E-01	-1.27E+00
4.00E+00	-2.30E-01	-1.51E+00	4.00E+00	-2.60E-01	-1.56E+00	4.00E+00	-2.45E-01	-1.48E+00	4.00E+00	-2.37E-01	-1.57E+00
4.50E+00	-2.22E-01	-1.86E+00	4.50E+00	-2.52E-01	-1.89E+00	4.50E+00	-2.37E-01	-1.85E+00	4.50E+00	-2.37E-01	-1.89E+00
5.00E+00	-2.15E-01	-2.26E+00	5.00E+00	-2.45E-01	-2.28E+00	5.00E+00	-2.37E-01	-2.26E+00	5.00E+00	-2.30E-01	-2.27E+00
5.50E+00	-2.15E-01	-2.72E+00	5.50E+00	-2.52E-01	-2.69E+00	5.50E+00	-2.37E-01	-2.72E+00	5.50E+00	-2.30E-01	-2.68E+00
6.00E+00	-2.08E-01	-3.22E+00	6.00E+00	-2.52E-01	-3.15E+00	6.00E+00	-2.37E-01	-3.23E+00	6.00E+00	-2.30E-01	-3.14E+00
Steel
t	x	y	t	x	y	t	x	y	t	x	y
0.00E+00	2.97E-01	-4.45E-02	0.00E+00	2.97E-01	-3.71E-02	0.00E+00	3.04E-01	-2.97E-02	0.00E+00	3.19E-01	-4.45E-02
5.00E-01	2.97E-01	-1.04E-01	5.00E-01	2.97E-01	-1.26E-01	5.00E-01	3.11E-01	-7.42E-02	5.00E-01	3.11E-01	-9.64E-02
1.00E+00	2.82E-01	-1.85E-01	1.00E+00	2.82E-01	-2.30E-01	1.00E+00	3.26E-01	-1.63E-01	1.00E+00	3.04E-01	-1.63E-01
1.50E+00	2.74E-01	-3.11E-01	1.50E+00	2.67E-01	-3.41E-01	1.50E+00	2.82E-01	-2.82E-01	1.50E+00	2.97E-01	-2.97E-01
2.00E+00	2.67E-01	-4.82E-01	2.00E+00	2.74E-01	-5.34E-01	2.00E+00	2.74E-01	-4.45E-01	2.00E+00	2.97E-01	-4.45E-01
2.50E+00	2.52E-01	-6.90E-01	2.50E+00	2.45E-01	-6.75E-01	2.50E+00	2.60E-01	-6.45E-01	2.50E+00	2.97E-01	-6.38E-01
3.00E+00	2.45E-01	-9.34E-01	3.00E+00	2.45E-01	-9.20E-01	3.00E+00	2.52E-01	-8.83E-01	3.00E+00	2.74E-01	-8.83E-01
3.50E+00	2.30E-01	-1.22E+00	3.50E+00	2.30E-01	-1.19E+00	3.50E+00	2.45E-01	-1.19E+00	3.50E+00	2.67E-01	-1.16E+00
4.00E+00	2.15E-01	-1.54E+00	4.00E+00	2.15E-01	-1.50E+00	4.00E+00	2.15E-01	-1.47E+00	4.00E+00	2.60E-01	-1.48E+00
4.50E+00	2.00E-01	-1.91E+00	4.50E+00	2.00E-01	-1.86E+00	4.50E+00	2.08E-01	-1.85E+00	4.50E+00	2.45E-01	-1.82E+00
5.00E+00	1.85E-01	-2.34E+00	5.00E+00	1.85E-01	-2.28E+00	5.00E+00	1.93E-01	-2.26E+00	5.00E+00	2.37E-01	-2.23E+00
5.50E+00	1.71E-01	-2.80E+00	5.50E+00	1.71E-01	-2.75E+00	5.50E+00	1.78E-01	-2.72E+00	5.50E+00	2.22E-01	-2.69E+00
6.00E+00	1.48E-01	-3.31E+00	6.00E+00	1.56E-01	-3.26E+00	6.00E+00	1.63E-01	-3.23E+00	6.00E+00	2.08E-01	-3.19E+00
"""

# ==================== Robust parser ====================
def parse_trials(lines):
    trials = [[] for _ in range(4)]
    for line in lines:
        tokens = line.replace('\t', ' ').split()
        if len(tokens) < 12: continue
        for i in range(4):
            try:
                t = float(tokens[3*i])
                y = float(tokens[3*i + 2])
                trials[i].append((t, y))
            except ValueError:
                continue
    processed = []
    for raw in trials:
        if len(raw) < 5: continue
        t_raw, y_raw = zip(*raw)
        t = np.array(t_raw) * TIME_SCALE
        y = np.array(y_raw)
        processed.append({'t': t, 'y': y})
    return processed

all_lines = [l for l in DATA.split('\n') if l.strip()]
mag_trials   = parse_trials(all_lines[2:all_lines.index('Steel')])
steel_trials = parse_trials(all_lines[all_lines.index('Steel')+2:])

# ==================== Averaging ====================
def average_trials(trials):
    all_t = np.sort(np.unique(np.concatenate([tr['t'] for tr in trials])))
    y_vals = np.array([np.interp(all_t, tr['t'], tr['y']) for tr in trials])
    return all_t, np.median(y_vals, axis=0), np.std(y_vals, axis=0, ddof=1)

tM, yM_mean, yM_std = average_trials(mag_trials)
tS, yS_mean, yS_std = average_trials(steel_trials)

# ==================== Nonlinear fit ====================
def model(t, v0, a):
    return v0 * t - 0.5 * a * t**2

poptM, pcovM = curve_fit(model, tM, yM_mean, p0=[-0.5, 9.8])
v0M, aM = poptM
perrM = np.sqrt(np.diag(pcovM))

poptS, pcovS = curve_fit(model, tS, yS_mean, p0=[-0.5, 9.8])
v0S, aS = poptS
perrS = np.sqrt(np.diag(pcovS))

# ==================== Plot with results in title ====================
plt.figure(figsize=(12, 7.5))

for tr in mag_trials:    plt.scatter(tr['t'], tr['y'], c='tab:blue', alpha=0.12, s=25)
for tr in steel_trials:  plt.scatter(tr['t'], tr['y'], c='tab:orange', alpha=0.12, s=25)

plt.fill_between(tM, yM_mean-yM_std, yM_mean+yM_std, color='tab:blue',  alpha=0.3, label='Magnet  mean ±σ')
plt.fill_between(tS, yS_mean-yS_std, yS_mean+yS_std, color='tab:orange', alpha=0.3, label='Steel   mean ±σ')

tt = np.linspace(0, 0.8, 600)
plt.plot(tt, model(tt, v0M, aM), color='navy',       lw=3.5, label=f'Magnet fit')
plt.plot(tt, model(tt, v0S, aS), color='darkorange', lw=3.5, label=f'Steel fit')

plt.axhline(0, color='k', lw=0.8)
plt.xlabel("Time (s) – 240→30 FPS corrected", fontsize=13)
plt.ylabel("Vertical position (m)", fontsize=13)
plt.grid(True, alpha=0.3)

title = (f"Magnet Down vs Steel Drop\n\n"
         f"Magnet: a = {aM:.3f} ± {perrM[1]:.3f} m/s²\n"
         f"Steel:  a = {aS:.3f} ± {perrS[1]:.3f} m/s²\n"
         f"Δa = {aS-aM:.3f} m/s²")
plt.title(title, fontsize=14, fontweight='bold', pad=20)

plt.legend(fontsize=12)
plt.tight_layout()
plt.show()
