"""Tracks the float in the night runs (side view, band y = 395 to 455 px, 10 fps).
Usage: python3 track_night.py <video> <name>   ->  <name>_track.csv: time (s), centre x, width, left, right (px)"""
import numpy as np, subprocess, sys, json
f,name=sys.argv[1],sys.argv[2]
W,H=1080,960; y0,y1=395,455; fps=10
p=subprocess.Popen(['ffmpeg','-v','error','-i',f,'-vf',f'fps={fps},crop={W}:{y1-y0}:0:{y0},format=gray','-f','rawvideo','-'],stdout=subprocess.PIPE)
out=[]; prev=None; i=0
while True:
    b=p.stdout.read(W*(y1-y0))
    if len(b)<W*(y1-y0): break
    g=np.frombuffer(b,np.uint8).reshape(y1-y0,W).astype(float)
    col=np.percentile(g,70,axis=0)            # bright if most of the band is bright
    col[:180]=0; col[1065:]=0
    thr=max(150,0.5*(np.percentile(col[180:1065],98)+np.percentile(col[180:1065],40)))
    m=col>thr
    # contiguous runs
    runs=[];s=None
    for x,v in enumerate(m):
        if v and s is None: s=x
        if (not v) and s is not None: runs.append((s,x-1)); s=None
    if s is not None: runs.append((s,len(m)-1))
    runs=[r for r in runs if r[1]-r[0]>60]
    if runs:
        if prev is None: r=max(runs,key=lambda r:r[1]-r[0])
        else: r=min(runs,key=lambda r:abs((r[0]+r[1])/2-prev))
        c=(r[0]+r[1])/2; prev=c; out.append((i/fps,c,r[1]-r[0],r[0],r[1]))
    else: out.append((i/fps,np.nan,np.nan,np.nan,np.nan))
    i+=1
a=np.array(out); np.savetxt(name+'_track.csv',a,delimiter=',',fmt='%.3f')
print(name,'frames',len(a),'start x %.0f w %.0f  end x %.0f w %.0f  nan %d'%(a[0,1],a[0,2],a[-1,1],a[-1,2],np.isnan(a[:,1]).sum()))
