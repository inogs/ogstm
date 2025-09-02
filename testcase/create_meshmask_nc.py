#! /usr/bin/python

# LOAD PACKAGES
import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

from interpolation43 import interpolate

def iseven(x):
    return not(bool(x&1))

# Script to create meshmask file

time = 1;
z_a  = 1;
y_a  = 1;
x_a  = 1;

#Earth Radius in meters
ER = 6371*10**3
PI = 3.141592653589793

def create_meshmask_nc(test):
    jpi=test['jpi']
    jpj=test['jpj']
    jpk=test['jpk']

    dx=test['dx']
    dy=test['dy']

#    np.float64 coastp(y, x) ;
    coastp = np.zeros((jpj,jpi),dtype=np.float64)


#    np.float64 fmask(time, z, y, x) ;
    fmask = np.ones((time,jpk,jpj,jpi),np.float64);
    fmask[0,:,1, 0]    = 0.;
    fmask[0,:,0, 1]    = 0.;
    fmask[0,:,-1,0]    = 0.;
    fmask[0,:,0,-1]    = 0.;

#    np.float64 gdept(time, z, y_a, x_a) ;
    filein             = 'KB'+'/gdept' + 'KB' + '.dat'
    gdeptTOT           = np.loadtxt(filein, dtype=np.float64);
    gdeptINT           = interpolate(gdeptTOT, jpk)
    gdept              = np.zeros((time,jpk,jpj,jpi),np.float64);
    for jj in range(jpj):
        for ji in range(jpi):
            gdept[0,0:jpk,jj,ji] = gdeptINT[0:jpk];

#    np.float64 gdepw(time, z, y_a, x_a) ;
    filein             = 'KB' + '/gdepw' + 'KB' + '.dat'
    gdepwTOT           = np.loadtxt(filein, dtype=np.float64);
    gdepwINT           = interpolate(gdepwTOT, jpk)
    gdepw              = np.zeros((time,jpk,y_a,x_a),np.float64);
    gdepw[0,0:jpk,0,0] = gdepwINT[0:jpk];


#    np.float64 glamt(time, z_a, y, x) ;
    glamt = np.ones((time,z_a,jpj,jpi),np.float64);
    if iseven(jpi):
        for ji in range(jpi):
            glamt[0,0,:,ji]=test['lon0']-jpi/2.*dx+dx/2. + dx*ji;
    else:
        for ji in range(jpi):
            glamt[0,0,:,ji]=test['lon0']-(jpi-1)/2.*dx + dx*ji;

#    np.float64 glamu(time, z_a, y, x) ;
    glamu = np.ones((time,z_a,jpj,jpi),np.float64);
    glamu = glamt+dx/2.;

#    np.float64 glamv(time, z_a, y, x) ;
    glamv = np.ones((time,z_a,jpj,jpi),np.float64);
    glamv = glamt;

#    np.float64 glamf(time, z_a, y, x) ;
    glamf = np.ones((time,z_a,jpj,jpi),np.float64);
    glamf = glamt+dx/2.;

#    np.float64 gphit(time, z_a, y, x) ;
    gphit = np.ones((time,z_a,jpj,jpi),np.float64);
    if iseven(jpj):
        for jj in range(jpj):
            gphit[0,0,jj,:]=test['lat0']-jpj/2.*dy+dy/2. + dy*jj;
    else:
        for jj in range(jpj):
            gphit[0,0,jj,:]=test['lat0']-(jpj-1)/2.*dy + dy*jj;

#    np.float64 gphiu(time, z_a, y, x) ;
    gphiu = np.ones((time,z_a,jpj,jpi),np.float64);
    gphiu = gphit

#    np.float64 gphiv(time, z_a, y, x) ;
    gphiv = np.ones((time,z_a,jpj,jpi),np.float64);
    gphiv = gphit + dy/2.

#    np.float64 gphif(time, z_a, y, x) ;
    gphif = np.ones((time,z_a,jpj,jpi),np.float64);
    gphif = gphit + dy/2.

#    np.float64 ff(time, z_a, y, x) ;
    ff = np.ones((time,z_a,jpj,jpi),np.float64);

#    float nav_lat(y, x) ;
    nav_lat = np.array(gphit,np.float32);

#    float nav_lev(z) ;
    nav_lev = np.array(gdept[0,0:jpk,0,0],np.float32);

#    float nav_lon(y, x) ;
    nav_lon = np.array(glamt,np.float32);

#    float time(time) ;
    timef = 1.;

#    short time_steps(time) ;
    time_steps = 1;

#    np.float64 e1f(time, z_a, y, x) ;
    e1f = np.ones((time,z_a,jpj,jpi),np.float64);
    for jj in range(jpj-1):
        for ji in range(jpi-1):
            dlam_di = glamf[0,0,jj,ji+1]-glamf[0,0,jj,ji]
            dphi_di = gphif[0,0,jj,ji+1]-gphif[0,0,jj,ji]
            e1f[0,0,jj,ji] = ER * PI/180. * np.sqrt( (dlam_di*np.cos(PI/180.*glamf[0,0,jj,ji]))**2 + dphi_di**2 )
    e1f[0,0,:,-1] =e1f[0,0,:,-2]
    e1f[0,0,-1,:] =e1f[0,0,-2,:]

#    np.float64 e1t(time, z_a, y, x) ;
    e1t = np.ones((time,z_a,jpj,jpi),np.float64);

    for jj in range(jpj-1):
        for ji in range(jpi-1):
            dlam_di = glamt[0,0,jj,ji+1]-glamt[0,0,jj,ji]
            dphi_di = gphit[0,0,jj,ji+1]-gphit[0,0,jj,ji]
            e1t[0,0,jj,ji] = ER * PI/180. * np.sqrt( (dlam_di*np.cos(PI/180.*glamt[0,0,jj,ji]))**2 + dphi_di**2 )
    e1t[0,0,:,-1] =e1t[0,0,:,-2]
    e1t[0,0,-1,:] =e1t[0,0,-2,:]

#    np.float64 e1u(time, z_a, y, x) ;
    e1u = np.ones((time,z_a,jpj,jpi),np.float64);

    for jj in range(jpj-1):
        for ji in range(jpi-1):
            dlam_di = glamu[0,0,jj,ji+1]-glamu[0,0,jj,ji]
            dphi_di = gphiu[0,0,jj,ji+1]-gphiu[0,0,jj,ji]
            e1u[0,0,jj,ji] = ER * PI/180. * np.sqrt( (dlam_di*np.cos(PI/180.*glamu[0,0,jj,ji]))**2 + dphi_di**2 )
    e1u[0,0,:,-1] =e1u[0,0,:,-2]
    e1u[0,0,-1,:] =e1u[0,0,-2,:]

#    np.float64 e1v(time, z_a, y, x) ;
    e1v = np.ones((time,z_a,jpj,jpi),np.float64);

    for jj in range(jpj-1):
        for ji in range(jpi-1):
            dlam_di = glamv[0,0,jj,ji+1]-glamv[0,0,jj,ji]
            dphi_di = gphiv[0,0,jj,ji+1]-gphiv[0,0,jj,ji]
            e1v[0,0,jj,ji] = ER * PI/180. * np.sqrt( (dlam_di*np.cos(PI/180.*glamv[0,0,jj,ji]))**2 + dphi_di**2 )
    e1v[0,0,:,-1] =e1v[0,0,:,-2]
    e1v[0,0,-1,:] =e1v[0,0,-2,:]

#    np.float64 e2f(time, z_a, y, x) ;
    e2f = np.ones((time,z_a,jpj,jpi),np.float64);

    for jj in range(jpj-1):
        for ji in range(jpi-1):
            dlam_dj = glamf[0,0,jj+1,ji]-glamf[0,0,jj,ji]
            dphi_dj = gphif[0,0,jj+1,ji]-gphif[0,0,jj,ji]
            e2f[0,0,jj,ji] = ER * PI/180. * np.sqrt( (dlam_dj*np.cos(PI/180.*glamf[0,0,jj,ji]))**2 + dphi_dj**2 )
    e2f[0,0,:,-1] =e2f[0,0,:,-2]
    e2f[0,0,-1,:] =e2f[0,0,-2,:]

#    np.float64 e2t(time, z_a, y, x) ;
    e2t = np.ones((time,z_a,jpj,jpi),np.float64)

    for jj in range(jpj-1):
        for ji in range(jpi-1):
            dlam_dj = glamt[0,0,jj+1,ji]-glamt[0,0,jj,ji]
            dphi_dj = gphit[0,0,jj+1,ji]-gphit[0,0,jj,ji]
            e2t[0,0,jj,ji] = ER * PI/180. * np.sqrt( (dlam_dj*np.cos(PI/180.*glamt[0,0,jj,ji]))**2 + dphi_dj**2 )
    e2t[0,0,:,-1] =e2t[0,0,:,-2]
    e2t[0,0,-1,:] =e2t[0,0,-2,:]

#    np.float64 e2u(time, z_a, y, x) ;
    e2u = np.ones((time,z_a,jpj,jpi),np.float64);

    for jj in range(jpj-1):
        for ji in range(jpi-1):
            dlam_dj = glamu[0,0,jj+1,ji]-glamu[0,0,jj,ji]
            dphi_dj = gphiu[0,0,jj+1,ji]-gphiu[0,0,jj,ji]
            e2u[0,0,jj,ji] = ER * PI/180. * np.sqrt( (dlam_dj*np.cos(PI/180.*glamu[0,0,jj,ji]))**2 + dphi_dj**2 )
    e2u[0,0,:,-1] =e2u[0,0,:,-2]
    e2u[0,0,-1,:] =e2u[0,0,-2,:]

#    np.float64 e2v(time, z_a, y, x) ;
    e2v = np.ones((time,z_a,jpj,jpi),np.float64);

    for jj in range(jpj-1):
        for ji in range(jpi-1):
            dlam_dj = glamv[0,0,jj+1,ji]-glamv[0,0,jj,ji]
            dphi_dj = gphiv[0,0,jj+1,ji]-gphiv[0,0,jj,ji]
            e2v[0,0,jj,ji] = ER * PI/180. * np.sqrt( (dlam_dj*np.cos(PI/180.*glamv[0,0,jj,ji]))**2 + dphi_dj**2 )
    e2v[0,0,:,-1] =e2v[0,0,:,-2]
    e2v[0,0,-1,:] =e2v[0,0,-2,:]


#   np.float64 e3t(time, z, y, x) ;
    filein             = 'KB'+'/e3t' + 'KB' + '.dat'
    e3tTOT             = np.loadtxt(filein, dtype=np.float64);
    e3tINT             = interpolate(e3tTOT, jpk)
    e3t                = np.zeros((time,jpk,jpj,jpi),np.float64);
    for jk in range(jpk):
        for jj in range(jpj):
            for ji in range(jpi):
                e3t[0,jk,jj,ji]   = e3tINT[jk];

#   np.float64 e3u(time, z, y, x) ;
    filein             = 'KB'+'/e3t' + 'KB' + '.dat'
    e3tTOT             = np.loadtxt(filein, dtype=np.float64);
    e3tINT             = interpolate(e3tTOT, jpk)
    e3u                = np.zeros((time,jpk,jpj,jpi),np.float64);
    for jk in range(jpk):
        for jj in range(jpj):
            for ji in range(jpi):
                e3u[0,jk,jj,ji]   = e3tINT[jk];

#   np.float64 e3v(time, z, y, x) ;
    filein             = 'KB'+'/e3t' + 'KB' + '.dat'
    e3tTOT             = np.loadtxt(filein, dtype=np.float64);
    e3tINT             = interpolate(e3tTOT, jpk)
    e3v                = np.zeros((time,jpk,jpj,jpi),np.float64);
    for jk in range(jpk):
        for jj in range(jpj):
            for ji in range(jpi):
                e3v[0,jk,jj,ji]   = e3tINT[jk];

#   np.float64 e3w(time, z, y_a, x_a) ;
    filein             = 'KB'+'/e3w' + 'KB' + '.dat'
    e3wTOT             = np.loadtxt(filein, dtype=np.float64);
    e3wINT             = interpolate(e3wTOT, jpk)
    e3w               = np.zeros((time,jpk,jpj,jpi),np.float64);
    for jk in range(jpk):
        for jj in range(jpj):
            for ji in range(jpi):
                e3w[0,jk,jj,ji]   = e3wINT[jk];

#    np.float64 tmask(time, z, y, x) ;
    tmask = np.ones((time,jpk,jpj,jpi),np.float64);
    tmask[0,:,0, :] =0.;
    tmask[0,:,:, 0] =0.;
    tmask[0,:,-1,:] =0.;
    tmask[0,:,:,-1] =0.;
    tmask[0,-1,:,:] =0.;

#    np.float64 umask(time, z, y, x) ;
    umask = np.ones((time,jpk,jpj,jpi),np.float64);
    umask[0,:,0, :] =0.;
    umask[0,:,:, 0] =0.;
    umask[0,:,-1,:] =0.;
    umask[0,:,:,-1] =0.;
    umask[0,:,:,-2] =0.;
    umask[0,-1,:,:] =0.;

#    np.float64 vmask(time, z, y, x) ;

    vmask = np.ones((time,jpk,jpj,jpi),np.float64);
    vmask[0,:,:, 0] =0.;
    vmask[0,:,0, :] =0.;
    vmask[0,:,-1,:] =0.;
    vmask[0,:,:,-1] =0.;
    vmask[0,:,-2,:] =0.;
    vmask[0,-1,:,:] =0.;

    ##############################################################
    # write meshmask netcdf file !
    ##############################################################
    os.system("mkdir -p " + test['Dir'].decode())

    outfile = test['Dir'].decode() + '/meshmask.nc';

    ncOUT=NC.netcdf_file(outfile,"w");

    ncOUT.createDimension('x',jpi);
    ncOUT.createDimension('y',jpj);
    ncOUT.createDimension('z',jpk);
    ncOUT.createDimension('time',time)

    ncOUT.createDimension('x_a',x_a);
    ncOUT.createDimension('y_a',y_a);
    ncOUT.createDimension('z_a',z_a);

    ncvar    = ncOUT.createVariable('coastp','d',('y','x'))                  ; ncvar[:] = coastp;
    ncvar    = ncOUT.createVariable('e1f'   ,'d',('time','z_a', 'y', 'x') )  ; ncvar[:] = e1f   ;
    ncvar    = ncOUT.createVariable('e1t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e1t   ;
    ncvar    = ncOUT.createVariable('e1u'   ,'d',('time','z_a', 'y', 'x')   ); ncvar[:] = e1u   ;
    ncvar    = ncOUT.createVariable('e1v'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e1v   ;
    ncvar    = ncOUT.createVariable('e2f'   ,'d',('time','z_a', 'y', 'x') )  ; ncvar[:] = e2f   ;
    ncvar    = ncOUT.createVariable('e2t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[:] = e2t   ;
    ncvar    = ncOUT.createVariable('e2u'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e2u   ;
    ncvar    = ncOUT.createVariable('e2v'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = e2v   ;
    ncvar    = ncOUT.createVariable('e3t_0' ,'d',('time','z'  , 'y', 'x'))   ; ncvar[:] = e3t   ;
    ncvar    = ncOUT.createVariable('e3u_0' ,'d',('time','z'  , 'y', 'x'))   ; ncvar[:] = e3u 
    ncvar    = ncOUT.createVariable('e3v_0' ,'d',('time','z'  , 'y', 'x'))   ; ncvar[:] = e3v
    ncvar    = ncOUT.createVariable('e3w_0' ,'d',('time','z'  , 'y', 'x'))   ; ncvar[:] = e3w 
    ncvar    = ncOUT.createVariable('ff'    ,'d',('time','z_a', 'y', 'x'))   ; ncvar[:] = ff    ;
    ncvar    = ncOUT.createVariable('fmask' ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = fmask ;
    ncvar    = ncOUT.createVariable('gdept' ,'d',('time','z', 'y', 'x'))     ; ncvar[:] = gdept ;
    ncvar    = ncOUT.createVariable('gdepw' ,'d',('time','z', 'y_a', 'x_a')) ; ncvar[:] = gdepw ;
    ncvar    = ncOUT.createVariable('glamf'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamf ;
    ncvar    = ncOUT.createVariable('glamt'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamt ;
    ncvar    = ncOUT.createVariable('glamu'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamu ;
    ncvar    = ncOUT.createVariable('glamv'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamv ;
    ncvar    = ncOUT.createVariable('gphif'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = glamf ;
    ncvar    = ncOUT.createVariable('gphit'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphit ;
    ncvar    = ncOUT.createVariable('gphiu'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphiu ;
    ncvar    = ncOUT.createVariable('gphiv'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[:] = gphiv ;
    ncvar    = ncOUT.createVariable('nav_lat','f',('y','x'))               ; ncvar[:] = nav_lat;
    ncvar    = ncOUT.createVariable('nav_lev' ,'f',('z',)) ; ncvar[:] = nav_lev;
    ncvar    = ncOUT.createVariable('nav_lon','f',('y','x'))                 ; ncvar[:] = nav_lon;
#	float time(time) ;
#	short time_steps(time) ;
    ncvar    = ncOUT.createVariable('tmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = tmask
    ncvar    = ncOUT.createVariable('umask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = umask
    ncvar    = ncOUT.createVariable('vmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = vmask
    ncOUT.close()

