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
#window_stride=points_t

tot_points_t=points_t+window_stride*(nwindows-1)
time_tot=Ts*tot_points_t

#%%

num_groups=1

num_users=get_num_participants()
t=np.arange(0,tot_points_t )*Ts


#groups=[ list( np.arange(1,num_users//2 +1 ) ) , list( np.arange(num_users//2 +1, num_users+1 ) )  ]

gsize_min=num_users//num_groups
gres=num_users%num_groups
assert gsize_min > 0 , "Too many groups for the number of participants"

usr_ids=range(1,num_users+1)
groups=[]
for g in range(num_groups):
    g1=g*gsize_min
    g2=g1+gsize_min+gres if g ==num_groups-1 else g1+gsize_min
    groups.append( usr_ids[g1:g2 ] )

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

#%%
Lx=dic_params["end_x"]
Ly=dic_params["end_y"]

#dic_pos={ 0: (0.0, 0.0) , 1: (Lx/2.0 , 0.0) , 2: ( Lx/2.0 , Ly/2.0 ) , 3: (Lx/2.0, Ly*3.0/4.0)  }
x0=Lx/2.0

#%% create points in space on the y axis
#we createthe points unevenly placed, logarithmically.
#generate log space
log_s=np.logspace(-1,1,nwindows-1)
log_s=np.append([0], log_s)
#normalize to the maximum number of the log space
log_s=log_s/np.max(log_s)
L_reg=Ly/2.0
#we scale our interval and invert it to make the center more highly sampled
scaled_log=L_reg+ log_s*L_reg
#invert it
ypos=scaled_log[:]


#%%
#ypos=np.floor(ypos)
df_space=pd.DataFrame({"x": x0, "y": ypos , "x_norm": x0/float(Lx) , 
                       "y_norm": ypos/float(Ly) })
df_space["win_ind"]=range( nwindows  )    
df_space.to_csv("space_points.csv")

#The order of the columns must be: 
# line_index, timepoint,  user, x,y,t

df_p1=pd.DataFrame({} )

for win in range(nwindows):
    ti=win*points_t
    tf=(win+1)*points_t
    tr=t[ti:tf ]
    df1=pd.DataFrame({} )
    #we add 1 because Fortran is 1 based: 
    df1["timep"]=np.arange(ti, tf).astype(int)+1
    df1["user"]=1
    df1["t"]=tr
    df1["x"]=x0
    df1["y"]=ypos[win]
    df_p1=df_p1.append(df1, sort=False, ignore_index=True)





#%%
df_all=df_p1.copy()
#raise Exception("Stopping")

#%%
for k in range(2, num_users+1):
    df2=df_p1.copy()
    df2["user"]=k
    df2["x"]=Lx/2.0
    df2["y"]=Ly/2.0
    df_all=pd.concat( [df_all, df2] ).reset_index(drop=True)

num_lines=len(df_all)

add_line="%d \t %d \t %d \n"%( num_lines , num_timepoints, num_users )
col_order=["timep","user", "x","y","t"]
#export to string in correct order
s=df_all[col_order].to_string()

#add the header line
file_string=add_line+s

#write gaze data to file
with open( file_out, "w" ) as fid:
    fid.write(file_string)

