#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:15:17 2020

@author: Yu Hong 
"""
from argparse import ArgumentParser
parser = ArgumentParser(description=__doc__)

parser.add_argument('imffolder', help='The name of the IMF Folder input file to read.',
                    type=str)
parser.add_argument('output_file', help='The name of the substorm output file to store.',
                    type=str)

# Get args from caller, collect arguments into a convenient object:
args = parser.parse_args()

# Import shit.  I needed a lot of shit this time.  
import numpy as np
import glob
from spacepy.pybats import gitm
from gitm_functs_sate import spheric_area, height_integ
from gitm_functs_sate import read_gitm, read_fac_gitm, sate_reader

#a = gitm.GitmBin('./dmsp/f16__001_t240510_082300.bin')
#%% define global parameters & constant
a=[]; b=[];c=[]; d=[]; e=[];f=[]; g=[]; h=[]
i=[]; j=[]; k=[]

#%% setting list & path for reading data
folder = args.imffolder

pathi="./"+folder+"/"
dirs = glob.glob(pathi+'*.bin')
    
# select one file for test or read files in order
for ifile,filename in enumerate(sorted(dirs)):  # dirs[::100]    
    if ifile <0:  
        continue
        
    print (ifile,filename)
        
# =============================================================================
# THIS part is used for polar region average of parameters of "HA"
# =============================================================================
    data = sate_reader(filename)
       
    lon = data[0]; lat = data[1]
    alt = data[2]/1000.; ne = data[4]; jh = data[5]
    rho = data[3]*10**12; # 9 for def 10 for fac 
    ef = data[3]
    vie = data[8]; vin = data[9]; viz = data[10]; 
    sgmp = data[13]; sgmh = data[14]
        
    hour = int(filename[25:27]); minu = int(filename[27:29])
    sec = int(filename[29:31])

    time_stampe = hour*3600 + minu*60 + sec
    #write this part into a function

    a.append(time_stampe); b.append(hour); c.append(minu)
    d.append(sec); e.append(lat); f.append(lon)
    g.append(alt); h.append(vie); i.append(vin); j.append(viz)
    data_file=open("./"+args.output_file+".txt",'w')
    #        data_file=open("./test_dmsp_read.txt",'w')
        
    for num,idata in enumerate(a):
        if num <0:
            continue
        data_file.write('%10d%5d%5d%6d%10.3f%10.3f%10.3f%15.5f%15.5f%15.5f' %(a[num],b[num],c[num],d[num],e[num],f[num],g[num],h[num],i[num],j[num])+'\n')
    data_file.close()





        
            
                           


