"""Tracks the float of the plastic-support runs (IMG_8765, IMG_8766) along the container.
The image is sampled at 4 frames per second; in a band at the waterline (rows 560-640) the mean brightness of each
column gives one line per frame. The static scene is removed by subtracting, column by column, the 5th percentile over
the record, and an edge of the float is followed from frame to frame, starting from a seed position.
Usage: python3 track_waterline.py <IMG_8765 video> <IMG_8766 video>   (writes tracks/IMG_876*_plastic_track.csv)"""
import numpy as np, subprocess, sys, os
from scipy.ndimage import gaussian_filter
here=os.path.dirname(os.path.abspath(__file__))
W,H,Y0,Y1=1080,960,560,640
def band(video):
    raw=subprocess.run(['ffmpeg','-v','error','-i',video,'-vf','fps=4,format=rgb24','-f','rawvideo','-'],capture_output=True).stdout
    n=len(raw)//(W*H*3); fr=np.frombuffer(raw,np.uint8).reshape(n,H,W,3)
    return np.stack([fr[k,Y0:Y1].astype(np.float32).mean(axis=(0,2)) for k in range(n)])
def follow(B,smooth,start,seeds,win):
    g=gaussian_filter(B,smooth); e=g-np.percentile(g,5,axis=0)
    a=gaussian_filter(e,(0.5,2.0)); d=np.empty_like(a)                    # slope along the line, central differences
    d[:,1:-1]=(a[:,2:]-a[:,:-2])/2; d[:,0]=a[:,1]-a[:,0]; d[:,-1]=a[:,-1]-a[:,-2]; n=len(e); out=[]
    for x0,sg in seeds:                        # sg = +1 rising edge, -1 falling edge
        xs=[None]*n
        def step(k,x):
            a,b=max(int(round(x))-win,1),min(int(round(x))+win+1,W-1); arr=sg*d[k]
            j=a+int(np.argmax(arr[a:b])); y1,y2,y3=arr[j-1],arr[j],arr[j+1]; den=y1-2*y2+y3
            return float(j+np.clip(0.5*(y1-y3)/den if den!=0 else 0,-1,1))
        x=x0
        for k in range(start,n): x=step(k,x); xs[k]=x
        x=xs[start]
        for k in range(start-1,-1,-1): x=step(k,x); xs[k]=x
        out.append(xs)
    return np.arange(n)/4.0, out
t,(xr,)=follow(band(sys.argv[1]),(1.0,2.5),0,[(668,+1)],10)              # one block: right edge of the float
np.savetxt(os.path.join(here,'tracks','IMG_8765_plastic_track.csv'),np.c_[t,xr],delimiter=',',fmt='%.4f',
           header='time (s), right edge of the float (px)')
t,(xa,xb)=follow(band(sys.argv[2]),(1.0,2.0),120,[(383,-1),(547,+1)],7)  # two blocks: two edges on the float, after its turn
np.savetxt(os.path.join(here,'tracks','IMG_8766_plastic_track.csv'),np.c_[t,xa,xb],delimiter=',',fmt='%.4f',
           header='time (s), left edge (px), right edge (px); the edges are followed from 30 s, after the turn')
print('tracks written')
