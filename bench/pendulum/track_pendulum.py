"""Tracks the pendulum bob by background subtraction (60 fps, 1080x1920 portrait).
Usage: python3 track_pendulum.py <video> <out.csv>   ->  time (s), x (px), y (px), blob size (px)"""
import numpy as np, subprocess, sys
from scipy import ndimage as ndi
f, out = sys.argv[1], sys.argv[2]
X0, Y0, CW, CH, S = 300, 560, 520, 640, 2      # crop in full-res pixels, downscale factor
w, h = CW//S, CH//S
raw = subprocess.run(['ffmpeg','-v','error','-i',f,'-vf',f'crop={CW}:{CH}:{X0}:{Y0},scale={w}:{h},format=gray','-f','rawvideo','-'],capture_output=True).stdout
fr = np.frombuffer(raw, np.uint8).reshape(-1, h, w)
fps = 60000/1001
bg = np.median(fr[::7].astype(np.float32), axis=0)
res = []
for i, a in enumerate(fr):
    d = np.abs(a.astype(np.float32) - bg)
    m = d > 28
    m = ndi.binary_opening(m, structure=np.ones((3,3)))     # drop the thin string
    lab, n = ndi.label(m)
    if n == 0: res.append((i/fps, np.nan, np.nan, 0)); continue
    sizes = ndi.sum(m, lab, range(1, n+1)); k = int(np.argmax(sizes)) + 1
    cy, cx = ndi.center_of_mass(d*(lab == k))
    res.append((i/fps, X0 + cx*S, Y0 + cy*S, sizes[k-1]))
np.savetxt(out, np.array(res), delimiter=',', fmt='%.4f')
r = np.array(res); print(out, len(r), 'found', np.isfinite(r[:,1]).mean(), 'median blob', np.nanmedian(r[:,3]))
