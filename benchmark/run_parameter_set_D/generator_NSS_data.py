#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 00:58:10 2020

@author: jorge
"""


import pandas as pd
import numpy as np


def get_num_participants():
    with open("num_participants.txt", "r") as fid:
        lines=fid.readline()
    line=lines.split()
    return int(line[0] )

def get_num_windows():
    with open("num_windows.txt", "r") as fid:
        lines=fid.readline()
    line=lines.split()
    return int(line[0] )

dic_params={ "samp_freq": 1.0, "points_t": 1, 
            "start_x": 1, "end_x": 1 ,
            "start_y": 1, "end_y": 1 ,
            "window_stride": 1 , 
            }

with open("parameters.txt","r") as fid:
    lines=fid.readlines()
lines_split=[ l.split() for l in lines if len(l.split() ) ]

pairs= [  (l[0],float( l[1] )  )   for  l in lines_split if dic_params.get(l[0],"")!=""    ] 

#Override defaults
for k,v in pairs:
    dic_params[k]=v


#%%

nwindows=get_num_windows()
fs=dic_params["samp_freq"]
Ts=1.0/fs
points_t=int( dic_params["points_t"] )
window_stride=int( dic_params["window_stride"] )

tot_points_t=points_t+window_stride*(nwindows-1)
time_tot=Ts*tot_points_t

#%%

num_groups=2

num_users=get_num_participants()
t=np.arange(0,tot_points_t )*Ts


groups=[ list( np.arange(1,num_users//2 +1 ) ) , list( np.arange(num_users//2 +1, num_users+1 ) )  ]
group_sizes=[ len(g) for g in groups   ]


#%%
with open("groups.txt","w") as fid_g:
    fid_g.write("%d\n"%num_groups)
    for g in range(num_groups):
        fid_g.write("%d\n"%group_sizes[g])
        users_str=" ".join([ "%d"%k for k in  groups[g] ] )        
        fid_g.write("%s\n"%users_str)

num_timepoints=len(t)
timepoints=range(1,num_timepoints+1)
file_out="gaze_data.dat"

x0=dic_params["start_x"]+dic_params["end_x"]*0.05
vx=0.8*(dic_params["end_x"]-x0)*1.0/(t[-1]-t[0])
x=x0+vx*t

y0=dic_params["end_y"]*0.01
vy=0
y=y0+vy*t

df=pd.DataFrame({} )

df["timep"]=timepoints
df["user"]=1
df["x"]=x
df["y"]=y
df["t"]=t
df_all=df.copy()
lane_width=(dic_params["end_y"]-y0)/(num_users+1)
#raise Exception("Stopping")

#%%
for k in range(2, num_users+1):
    df2=df.copy()
    df2["user"]=k
    df2["y"]=df2["y"]+lane_width*(k-1)
    df2["x"]=df2["x"]+0*(k-1)
    df_all=pd.concat( [df_all, df2] ).reset_index(drop=True)

num_lines=len(df_all)

add_line="%d \t %d \t %d \n"%( num_lines , num_timepoints, num_users )
s=df_all.to_string()

file_string=add_line+s

with open( file_out, "w" ) as fid:
    fid.write(file_string)

