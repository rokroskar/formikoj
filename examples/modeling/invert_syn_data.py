import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import TravelTimeManager

# read the pck file
srtData = pg.DataContainer("./out/syn_data/syn_data_tt.pck", 
                           sensorTokens="s g",
                           removeInvalid=True)
abserr = 0.0035 # absolute error of traveltimes (s)
srtData.set('err', np.ones(srtData.size()) * abserr)

# create TravelTimeManger object
mgr = TravelTimeManager(srtData)

# create mesh
plc = mt.createParaMeshPLC(srtData.sensorPositions(),
                           paraDX=0.25, # number of nodes between stations (sensor positions)
                           paraDepth=15, # model depth (m) (0 --> estimated)
                           paraBoundary=2, # margin for parameter domain 
                                           # (sensor spacing)
                           paraMaxCellSize=2, # maximum cell area (m^2)
                           boundary=0) # boundary around the parameter domain
                                       # no need to worry about this here
mesh = mt.createMesh(plc,
                     quality=32, # 2D triangle quality
                                 # recommended values: 30 - 34
                     smooth=(1, 10)) # smoothing alogrithm - no need to change this

# set the mesh
mgr.setMesh(mesh,
            secNodes=3) # secondary nodes to refine the mesh to
                        # permit a more accurate ray tracing

# define the inversion parameters
lam = 50 # use different values
zWeight = .5 # isotropic smoothing - use for all inversions
                        
# run inversion
vest = mgr.invert(vTop=500, # velocity at the top of the initial model (m/s)
                  vBottom=5000, # velocity at the bottom of the initial model (m/s)
                  lam=lam, # regularization parameter
                  zWeight=zWeight, # 1 - same smoothing in x and z direction
                                   # 0.25 - four times stronger smoothing in 
                                   #        in horizontal direction
                  robustData=False, verbose=True) # use iterative re-weighting

# Consider the ray coverage in the plot, i.e.,
# blank areas in the model that should not be
# interpreted.
scov = mgr.standardizedCoverage()
vest = np.ma.masked_where(scov == 0, vest)

# plot the inversion result
cMin = 750
cMax = 4000

fig, ax = pg.plt.subplots(1, figsize=(8, 5))

_, cb = pg.show(mesh, # mesh
                vest, # resolved velocities
                coverage=scov, # blank areas with low coverage
                # ~ cMap='inferno_r', # colormap
                cMap='Spectral', # colormap
                # ~ cMap='plasma', # colormap
                logScale=False, # illustrate log-scaled velocities
                cMin=cMin, # lower limit of colormap
                cMax=cMax, # upper limit of colormap
                colorBar=True, # show colorbar
                orientation='horizontal', # orientation of the colorbar
                ax=ax, # axis object to plot in
                fitView=True, # automatic aspect ratio
                showBoundary=False) # plot mesh boundaries
                
# format plot
ax.set_xlabel('x (m)')
ax.set_ylabel('z (m)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
cb.set_label('v$_p$ ($\mathregular{m\/s^{-1}}$)', fontsize=10)

pg.plt.show()
