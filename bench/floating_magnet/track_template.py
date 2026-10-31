"""Tracks the float in the basement runs by edge-template matching (10 fps; template = float in the first frame).
Usage: python3 track_template.py <video> <out.csv> <template x0> <template x1>   (IMG_8751: 375 690; IMG_8750: 235 575)
 -> time (s), centre x (px), match score; IMG_8750 matches poorly (score below 0.5) and is described from frames in the paper"""
import numpy as np, subprocess, sys
from scipy.ndimage import sobel, gaussian_filter
f, out, tx0, tx1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
W, Y0, Y1, FPS = 1080, 575, 745, 10
raw = subprocess.run(['ffmpeg','-v','error','-i',f,'-vf',f'fps={FPS},crop={W}:{Y1-Y0}:0:{Y0},format=gray','-f','rawvideo','-'],capture_output=True).stdout
fr = np.frombuffer(raw, np.uint8).reshape(-1, Y1-Y0, W).astype(float)
def feat(b):
    g = gaussian_filter(b, 2); e = np.abs(sobel(g, axis=1)) + np.abs(sobel(g, axis=0)); return e
tmpl = feat(fr[0])[:, tx0:tx1]; tw = tx1 - tx0
tz = (tmpl - tmpl.mean())/tmpl.std()
res = []; prev = (tx0+tx1)/2
for i, b in enumerate(fr):
    e = feat(b); best = (-9, 0)
    for x in range(0, W - tw, 3):
        p = e[:, x:x+tw]; s = p.std()
        if s == 0: continue
        c = ((p - p.mean())/s * tz).mean()
        if c > best[0]: best = (c, x)
    res.append((i/FPS, best[1] + tw/2, best[0]))
np.savetxt(out, np.array(res), delimiter=',', fmt='%.3f'); print(out, len(res))
