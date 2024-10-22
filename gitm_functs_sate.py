#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 17:36:20 2020

@author: yxh5920
"""
# =============================================================================
# THESE are functions for processing data from GITM
# =============================================================================
import numpy as np
from spacepy.pybats import gitm
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#from funct_plot import make_mid
 
def read_gitm(filename):
    """
    this is used for reading GITM 3DALL files
    """
#    mdir = '/Users/yxh5920/Desktop'
#    a = gitm.GitmBin(mdir + filename)
    a = gitm.GitmBin(filename)
    p=180.0/np.pi
    #lon = a['Longitude'][2:-2,14:62,13]*p
    lon = a['Longitude'][2:-2,2:-2,13]*p
#    lat = p*a['Latitude'][2:-2,14:62,13]
    lat = p*a['Latitude'][2:-2,2:-2,13]
    par = a['Potential'][2:-2,2:-2,13]
    ae = a['ElectronAverageEnergy'][2:-2,2:-2,13]
    ef = a['Energyflux'][2:-2,2:22,13] # 2:52 when use step1+step2
    alt = a['Altitude'][5,5,2:52] # 2:52  13
#    alt = a['Altitude'][5,5,2:72] # 2:52  13
    joule = a['Joule heating'][2:-2,2:-2,2:52] # 13 for E-region 27 for F region JH only E-egion 36 for 400km
    jh = a['Joule heating'][2:-2,54:-2,2:52] # south-hemisphere 2:22 north: 54:-2 full 2:-2
    jh_s = a['Joule heating'][2:-2,2:22,2:52]
#    Ne = a['e-'][2:-2,2:-2,13]
    Ne = a['e-'][2:-2,2:-2,2:52]#
#    Ne = a['e-'][2:-2,2:-2,2:72]# 
    Rho = a['Rho'][2:-2,2:-2,36]# 13 - 123km 27 -254km 30 -300km 36 -400km 40-465km
    Ve = a['V!Dn!N (east)'][2:-2,2:-2,13]
    Vn = a['V!Dn!N (north)'][2:-2,2:-2,13] # 23 for 200km
    Vz = a['V!Dn!N (up)'][2:-2,2:-2,13] 
    # Ve = a['V!Dn!N (east)'][2:-2,14:62,35]
    # Vn = a['V!Dn!N (north)'][2:-2,14:62,35] # 23 for 200km
    # Vz = a['V!Dn!N (up)'][2:-2,14:62,35] 
    # Vin = a['V!Di!N (north)'][2:-2,14:62,35]; Vie = a['V!Di!N (east)'][2:-2,14:62,35]
    # Viz = a['V!Di!N (up)'][2:-2,14:62,35];
    Vin = a['V!Di!N (north)'][2:-2,2:-2,13]; Vie = a['V!Di!N (east)'][2:-2,2:-2,13]
    Viz = a['V!Di!N (up)'][2:-2,2:-2,13];
    # gedy_pot = a['Gedy_pot'][2:-2,2:-2,13];
    
    #Efield = a['EField'][2:-2,2:-2,13];
    #sgm_P = a['Pedersen Conductance'][2:-2,2:-2,13]
    # sgm_H = a['Hall Conductance'][2:-2,54:-2,13]
    
    # convect dmarray into array (dmarray = datamodel)
    x = np.array(lon.tolist())
    y = np.array(lat.tolist())
    z = np.array(par.tolist())
    w = np.array(ef.tolist())
    u = np.array(ae.tolist())
    v = np.array(joule.tolist())
    h = np.array(alt.tolist())
    j = np.array(jh.tolist())
    ne = np.array(Ne.tolist())
    rho = np.array(Rho.tolist())
    ve = np.array(Ve.tolist())
    vn = np.array(Vn.tolist())
    vz = np.array(Vz.tolist())
    vin = np.array(Vin.tolist()); vie = np.array(Vie.tolist())
    viz = np.array(Viz.tolist())
#    sjh = np.array(jh_s.tolist())
   # gedy_pot = np.array(gedy_pot.tolist())
    #efield = np.array(Efield.tolist())
    #sgmp = np.array(sgm_P.tolist())
    # sgmh = np.array(sgm_H.tolist())
    
    return x,y,z,w,u,v,h,j,ne,rho,ve,vn,vz,vie,vin,viz#,gedy_pot#,sgmp#,efield #, #,sjh#, #,,sgmh

def height_integ(j,h):
    """
    This is used for calculating the height-integral of Joule Heating
    """
    nLat = j.shape[1]; nLon = j.shape[0]
    int_value_2 = np.zeros((nLon,nLat)) # *np.nan
    del_ht = (h[1:]-h[:-1])
    ave_jh = (j[:,:,1:]+j[:,:,:-1])/2.0
    
    for ii in range(nLat): 
        for jj in range(nLon):
            int_value_2[0:jj+1,0:ii+1] = np.dot((ave_jh[0:jj+1,0:ii+1,:]),del_ht)
#            int_value_add_2 = (int_value_2[0]+int_value_2[71])/2
#            int_value_new_2 = np.vstack((int_value_add_2,int_value_2,int_value_add_2)) # this is the heighted-integrated JH from 41.25 to 88.75 degree
            
    return int_value_2


def spheric_integ(nLat,nLon,delt_lat,int_value_2):
    """
     THIS is used for calculating the hemispheric-integral of Joule Heating
    """
    band_area = np.zeros((nLat,1)) # from lat = 88.75 to 41.25
    for ilat in range(nLat):
        band_area[ilat] = np.sin((88.75-delt_lat*ilat)*np.pi/180)-np.sin((88.75-delt_lat*(ilat+1))*np.pi/180)
        jh_theta = np.sum(int_value_2,axis=0)/nLon # sum of height-integrated JH in each latitude band 
        jh_theta2 = jh_theta[0:20] # select region: 45-90 from 41-90/ 20 for default 10 for new
        jh_sphere = ((np.dot(jh_theta2[::-1],band_area[0:20])))*(2*(np.pi)*(110000+6371000)**2)/(10**9)
        # if for southern hemisphere: delect [::-1]!!!!!
        
    return jh_sphere


def spheric_area(nLat,nLon,delt_lat,int_value_2):
    """
     THIS is used for calculating the hemispheric-area integral of Joule Heating
    """
    length = np.zeros((nLat,1)) 
    band_area2 = np.zeros((nLat,1)) 
    for iilat in range(nLat):
        jh_theta = np.sum(int_value_2,axis=0)/nLon
        jh_theta2 = jh_theta[2:20] # select region: 45-90 from 41-90
        length[iilat] = ((2*(np.pi)*(110000+6371000)))*(np.sin((iilat*2.5)*np.pi/180))
        height = (110000+6371000)*(np.pi/180)*delt_lat
        band_area2[iilat] = length[iilat]*height
        jh_sphere2 = ((np.dot(jh_theta2[::-1], band_area2[1:-1])))/(10**9)
        
    return jh_sphere2

# =============================================================================
#     delt_lon = 5; delt_lat = 2.5;
# =============================================================================

def horizon_wind(m,n,v_e,v_n):
    """
    THIS is used for calculating the magnitude of the horizontal wind
    v_e --> Zonal wind; v_n --> Meridional wind
    """
    v_t = np.zeros((m,n))
    
    for ilon in range(len(v_e[:])):
        for ilat in range(len(v_e[1])):
            v_t[ilon][ilat] = np.sqrt(np.square(v_n[ilon][ilat])+np.square(v_e[ilon][ilat]))
            
    return  v_t    

def cat2pol(m,n,theta,v_e,v_n):
    """
    THIS is used for projecting the Cartesian coordinate to Polar coordinate
    v_e --> Zonal wind; v_n --> Meridional wind
    """
    v_e_pl = np.zeros((m,n)); v_n_pl = np.zeros((m,n))
    angle = theta[:,0]
    for i in range(len(angle)):
        if i<(len(angle)/2): 
            
            v_e_pl[i][:] = (v_e[i][:])*(np.cos(angle[i]))-(v_n[i][:])*(np.sin(angle[i])) # -0.5*np.pi
            v_n_pl[i][:] = (v_n[i][:])*(np.cos(angle[i]))+(v_e[i][:])*(np.sin(angle[i])) # + --> -
        
        else:
            v_e_pl[i][:] = (v_e[i][:])*(np.cos(angle[i]))-(v_n[i][:])*(np.sin(angle[i]))
            v_n_pl[i][:] = (v_n[i][:])*(np.cos(angle[i]))+(v_e[i][:])*(np.sin(angle[i])) # + --> -
    
    return v_e_pl, v_n_pl
    

# =============================================================================
# This is used for GITM.bin data rerange from ghost cells into full cells
# =============================================================================
def bin2ns(lon,lat,parm):
    """
    This is used for separate GITM.bin  data into N-S hemispheres
    """
    p = 180.0/np.pi
    lon_add_1 = (abs(lon[71])*0)
    lon_add_2 = (abs(lon[71])*2+5)/2
    lon_new = np.vstack((lon_add_1,lon,lon_add_2))
    
    parm_add = (parm[0]+parm[71])/2
    parm_new = np.vstack((parm_add,parm,parm_add))
    
    lat_add_1 = (lat[0])
    lat_new = np.vstack((lat_add_1,lat,lat_add_1))
    
    lon0 = lon_new/p; theta = lon0[:,52:72]; radiu = lat_new[:,52:72]
    
    nparm = parm_new[:,52:72]; 
    parm0 = parm_new[:,::-1]; sparm = parm0[:,52:72]
    
    return theta, radiu, nparm, sparm


# =============================================================================
# This is used for reading gitm_data altitudinal profile
# =============================================================================
def read_gitm_alt(filename):
    """
    this is used for reading GITM 3DALL files
    """
    a = gitm.GitmBin(filename)
    p=180.0/np.pi
    #define the parameters       
    # lon = a['Longitude'][2:-2,2:-2,13]*p
    # lat = p*a['Latitude'][2:-2,2:-2,13]
    lon = a['m-lon'][2:-2,2:-2,13]*p
    lat = p*a['m-lat'][2:-2,2:-2,13]
    alt = a['Altitude'][2:-2,2:-2,2:52]
    joule = a['Joule Heating'][2:-2,2:-2,2:52] 
    jh = a['Joule Heating'][2:-2,2:-2,2:52] 
    Ne = a['e-'][2:-2,2:-2,2:52]
    Rho = a['Rho'][2:-2,2:-2,30:52]
    Ve = a['V!Dn!N (east)'][2:-2,2:-2,2:52]
    Vn = a['V!Dn!N (north)'][2:-2,2:-2,2:52] 
    Vz = a['V!Dn!N (up)'][2:-2,2:-2,2:52];
    Vz = a['V!Dn!N (up)'][2:-2,2:-2,2:52] 
    Vin = a['V!Di!N (north)'][2:-2,2:-2,2:52]
    Vie = a['V!Di!N (east)'][2:-2,2:-2,2:52]
    Viz = a['V!Di!N (up)'][2:-2,2:-2,2:52];
    Efield = a['EField'][2:-2,2:-2,2:52];
    sgm_P = a['Pedersen Conductance'][2:-2,2:-2,2:52]
    sgm_H = a['Hall Conductance'][2:-2,2:-2,2:52]
    
    for var in [lon,lat,alt,joule,jh,Ne,Rho,Ve,Vn,Vz,Vin,Vie,Viz,Efield,sgm_P,sgm_H]:
        var = np.array(var.tolist())
    
    return lon,lat,alt,joule,jh,Ne,Rho,Ve,Vn,Vz,Vin,Vie,Viz,Efield,sgm_P,sgm_H

# =============================================================================
# This is used for reading FAC-module GITM simulated data
# =============================================================================
def read_fac_gitm(filename):
    """
    this is used for reading GITM 3DALL files
    """
    a = gitm.GitmBin(filename)
    p=180.0/np.pi
    # define the parameters       
    lon = a['Longitude'][2:-2,2:-2,13]*p
    lat = p*a['Latitude'][2:-2,2:-2,13]
    pot = a['Potential'][2:-2,2:-2,13]
    ae = a['ElectronAverageEnergy'][2:-2,2:-2,13]
    ef = a['Energyflux'][2:-2,54:-2,13] # 2:52 when use step1+step2
    alt = a['Altitude'][5,5,2:52] # 2:52  13 old-grid: 2:52/new-grid: 2:72
    jh = a['Joule heating'][2:-2,2:-2,2:52] # 13 for E-region 27 for F region JH only E-egion 36 for 400km
    jh_n = a['Joule heating'][2:-2,54:-2,2:52] # south-hemisphere 2:22 north: 54:-2 full 2:-2 28:-2
    jh_s = a['Joule heating'][2:-2,2:22,2:52] # old-grid: 54, 22; -->12
    Ne = a['e-'][2:-2,2:-2,13]
    Rho = a['Rho'][2:-2,2:-2,30]# 13 - 123km 27 -254km 30 -300km 36 -400km 40-465km
    Ve = a['V!Dn!N (east)'][2:-2,2:-2,30]
    Vn = a['V!Dn!N (north)'][2:-2,2:-2,30] # 23 for 200km
    Vz = a['V!Dn!N (up)'][2:-2,2:-2,13] 
    Vin = a['V!Di!N (north)'][2:-2,2:-2,13]; Vie = a['V!Di!N (east)'][2:-2,2:-2,13]
    Viz = a['V!Di!N (up)'][2:-2,2:-2,13];
    
    # SH: 2:12 NH: 28:-2
    gedy_fac = a['Gedy_fac'][2:-2,28:-2,13]; gedy_fac1 = a['Gedy_fac1'][2:-2,2:-2,13]
    gedy_pot = a['Gedy_pot'][2:-2,2:-2,13]; gedy_pot1 = a['Gedy_pot1'][2:-2,2:-2,13]
    gedy_dfac = a['Gedy_dfac'][2:-2,2:-2,13]; 
    
    Efield = a['E.F. Magnitude'][2:-2,2:-2,13]; Ef_n = a['E.F. North'][2:-2,2:-2,13]
    Ef_e = a['E.F. East'][2:-2,2:-2,13]; Ef_z = a['E.F. Vertical'][2:-2,2:-2,13]
    sgm_P = a['PedersenConductance'][2:-2,2:-2,13]; sgm_H = a['HallConductance'][2:-2,2:-2,13]
    
    # convect dmarray into array (dmarray = datamodel)
    lon = np.array(lon.tolist())
    lat = np.array(lat.tolist())
    pot = np.array(pot.tolist())
    ef = np.array(ef.tolist())
    ae = np.array(ae.tolist())
    jh = np.array(jh.tolist()); jh_n = np.array(jh_n.tolist()); jh_s = np.array(jh_s.tolist())
    alt = np.array(alt.tolist())
    ne = np.array(Ne.tolist())
    rho = np.array(Rho.tolist())
    ve = np.array(Ve.tolist()); vn = np.array(Vn.tolist()); vz = np.array(Vz.tolist())
    vin = np.array(Vin.tolist()); vie = np.array(Vie.tolist()); viz = np.array(Viz.tolist())
    gedy_pot = np.array(gedy_pot.tolist()); gedy_pot1 = np.array(gedy_pot1.tolist())
    gedy_fac = np.array(gedy_fac.tolist()); gedy_fac1 = np.array(gedy_fac1.tolist())
    gedy_dfac = np.array(gedy_dfac.tolist()); efield = np.array(Efield.tolist())
    ef_n = np.array(Ef_n.tolist()); ef_e = np.array(Ef_e.tolist()); ef_z = np.array(Ef_z.tolist())
    sgmp = np.array(sgm_P.tolist()); sgmh = np.array(sgm_H.tolist())
    
    return lon,lat,pot,ef,ae,jh,jh_n,jh_s,alt,ne,rho,ve,vn,vz,vie,vin,viz,gedy_pot,\
    gedy_pot1,gedy_fac,gedy_fac1,gedy_dfac,efield,ef_n,ef_e,ef_z,sgmp,sgmh


# =============================================================================    
def sate_reader(filename):
    """
    this is used for reading GITM #Satellite files
    """
    a = gitm.GitmBin(filename)
    p=180.0/np.pi
    #%% define the parameters
    lon = a['Longitude'][13]*p; lat = p*a['Latitude'][13]; 
    alt = a['Altitude'][34]
    rho = a['Rho'][36]; ne = a['e-'][36]
    ve = a['V!Dn!N (east)'][36]; vn = a['V!Dn!N (north)'][36] 
    vz = a['V!Dn!N (up)'][36]; vin = a['V!Di!N (north)'][34]
    vie = a['V!Di!N (east)'][34]; viz = a['V!Di!N (up)'][34]
    ae = a['ElectronAverageEnergy'][36];# pot = a['Potential'][:,:,13]
    ae = a['ElectronAverageEnergy'][36]; ef = a['Energyflux'][36]
    #jh = a['Joule heating'][:,:,13]; Efield = a['E.F. Magnitude'][:,:,13]
    #sgm_P = a['PedersenConductance'][:,:,13]; sgm_H = a['HallConductance'][:,:,13]
    #Ef_n = a['E.F. North'][2:-2,2:-2,13]; Ef_e = a['E.F. East'][:,:,13]
    #Ef_z = a['E.F. Vertical'][:,:,13]
    # convect dmarray into array (dmarray = datamodel)
#    for var in [lon,lat,alt,rho,ne,ve,vn,vz,vie,vin,viz,ae,pot,ae,ef,jh,Efield,sgm_P, \
#                sgm_H,Ef_e,Ef_n,Ef_z]:

    for var in [lon,lat,alt,rho,ne,ve,vn,vz,vie,vin,viz,ae,ae,ae,ef]:
        var = np.array(var.tolist())
    return lon,lat,alt,rho,ne,ve,vn,vz,vie,vin,viz,ae,ae,ae,ef#,jh,Efield,sgm_P, \
        #sgm_H,Ef_e,Ef_n,Ef_z # ae pot ae

# ======================================================================
# This is used for reading 2DUSR* GITM simulated files
# ======================================================================    
def read_2dusr(filename):
    a = gitm.GitmBin(filename)
    p=180.0/np.pi
    pot = a['Potential'][2:-2,2:-2]
    ef = a['EFlux'][2:-2,2:-2]
    jh = a['intJH'][2:-2,2:22] # 2:52 when use step1+step2  
    ef_n = a['EFlux'][2:-2,54:-2]; ef_s = a['EFlux'][2:-2,2:22]

    pot = np.array(pot.tolist())
    ef = np.array(ef.tolist()); ef_n = np.array(ef.tolist())
    jh = np.array(jh.tolist()); ef_s = np.array(ef.tolist())

    return pot,ef,jh,ef_n,ef_s









    
