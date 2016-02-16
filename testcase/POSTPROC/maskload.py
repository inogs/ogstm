import scipy.io.netcdf as NC
import numpy as np
import os, sys

def getDepthIndex(nav_lev, lev_mask):
    jk_m = 0
    for jk in range(nav_lev.__len__()):
        if nav_lev[jk] <= lev_mask:
            jk_m = jk
    return jk_m



maskfile    = "../TEST02/wrkdir/MODEL/"  + "meshmask.nc"
#maskfile    = os.getenv("MODDIR") + "meshmask.nc"

M=NC.netcdf_file(maskfile,"r")

jpi=M.dimensions['x']
jpj=M.dimensions['y']
jpk=M.dimensions['z']

tmask      = (M.variables['tmask'].data[0,:,:,:]).astype(np.bool).copy()

nav_lev    =  np.abs(M.variables['nav_lev'].data.copy())
e3t        =  np.abs(M.variables['e3t'].data.copy())

Lon        =  M.variables['nav_lon'].data[:,:].copy()
Lat        =  M.variables['nav_lat'].data[:,:].copy()

tk_m       = getDepthIndex(nav_lev,200.0)

mask200_2D = tmask[tk_m,:,:].copy()
mask200_3D = np.zeros((jpk,jpj,jpi),dtype=np.bool)

for i in range(jpk):
    mask200_3D[i,:,:]=mask200_2D

COASTNESS_LIST=['coast','open_sea','everywhere']
struct=[]
for coast in COASTNESS_LIST:
    struct.append((coast,np.bool))
    
COASTNESS = np.ones((jpk,jpj,jpi),dtype=struct)
COASTNESS['coast']     = ~mask200_3D;
COASTNESS['open_sea']  =  mask200_3D;

DEPTHlist      =['shallow','deep']
struct=[]
for depth in DEPTHlist:
    struct.append((depth,np.bool))

DEPTH  = np.zeros((jpk,jpj,jpi),dtype=struct)
tk_1   = getDepthIndex(nav_lev,200.0)+1

bathy  = nav_lev[tmask.sum(0)]
