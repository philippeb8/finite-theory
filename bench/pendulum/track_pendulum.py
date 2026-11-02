"""Tracks the pendulum bob by background subtraction (60 fps, 1080x1920 portrait).
The background is the median of the frames within +-3 s, renewed every second, so that it follows the small shift
of the camera (about 2 pixels) between 52.6 and 54.5 s of IMG_8754; a single background lost the ball after 54 s.
Usage: python3 track_pendulum.py <video> <out.csv>   ->  time (s), x (px), y (px), blob size (px)"""
import numpy as np, subprocess, sys
from scipy import ndimage as ndi
f, out = sys.argv[1], sys.argv[2]
X0, Y0, CW, CH, S = 300, 560, 520, 640, 2      # crop in full-res pixels, downscale factor
w, h = CW//S, CH//S
raw = subprocess.run(['ffmpeg','-v','error','-i',f,'-vf',f'crop={CW}:{CH}:{X0}:{Y0},scale={w}:{h},format=gray','-f','rawvideo','-'],capture_output=True).stdout
fr = np.frombuffer(raw, np.uint8).reshape(-1, h, w)
fps = 60000/1001; n = len(fr); B = int(round(fps)); half = 3*B
bgs = {}
def background(i):                                 # median of the frames within +-3 s of the second holding frame i
    k = i//B
    if k not in bgs:
        c = k*B + B//2; bgs[k] = np.median(fr[max(0, c-half):min(n, c+half):4].astype(np.float32), axis=0)
    return bgs[k]
res = []
for i, a in enumerate(fr):
    d = np.abs(a.astype(np.float32) - background(i))
    m = d > 28
    m = ndi.binary_opening(m, structure=np.ones((3,3)))     # drop the thin string
    lab, nl = ndi.label(m)
    if nl == 0: res.append((i/fps, np.nan, np.nan, 0)); continue
    sizes = ndi.sum(m, lab, range(1, nl+1)); k = int(np.argmax(sizes)) + 1
    cy, cx = ndi.center_of_mass(d*(lab == k))
    res.append((i/fps, X0 + cx*S, Y0 + cy*S, sizes[k-1]))
np.savetxt(out, np.array(res), delimiter=',', fmt='%.4f')
r = np.array(res); print(out, len(r), 'found', np.isfinite(r[:,1]).mean())
