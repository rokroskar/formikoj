import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import TravelTimeManager

from formikoj import SeismicWaveformModeler

#### CREATE MODEL ####

# create the model geometry

top = mt.createPolygon([[0, 0], [94, 0], 
                        [94, -5.5], [72, -5.5], 
                        [20, -2], [0, -2]],
                       isClosed=True, marker=1, area=0.1)

infill = mt.createPolygon([[20, -2], [72, -5.5], 
                           [70, -11], [22, -11]],
                          isClosed=True, marker=2, area=0.1)

bottom = mt.createPolygon([[0, -2], [20, -2], 
                           [22, -11], [70, -11], 
                           [72, -5.5], [94, -5.5], 
                           [94, -15], [0, -15]],
                          isClosed=True, marker=3, area=0.1)

geom = top + infill + bottom

# ~ pg.show(geom)
# ~ plt.show()

# add shot and receiver stations
nstats = 48
spc = 2

stations = np.vstack((np.linspace(0, (nstats-1)*spc, nstats), 
                      np.zeros(nstats))).T

for p in stations:
    geom.createNode(p)

# create mesh for the finite element modeling
mesh = mt.createMesh(geom, quality=34)
mesh.save('./in/mesh.bms')

# assign velocity values to the different model regions
velmap = [[1, 750],
          [2, 2500],
          [3, 4000]]
          
# ~ # show the final parameterized model
# ~ pg.show(mesh, data=velmap, label=pg.unit('vel'), showMesh=True)
# ~ plt.show()

#### CREATE SYNTHETIC TRAVEL TIME DATA ####

swm = SeismicWaveformModeler('.')

swm.load('mesh')
swm.load('config')
# ~ swm.create('waveforms')
swm.create('traveltimes')

# ~ #### INVERT SYNTHETIC TRAVEL TIMES ####

# ~ # read the pck file
# ~ srtData = pg.DataContainer("./out/syn_data/syn_data_tt.pck", 
                           # ~ sensorTokens="s g",
                           # ~ removeInvalid=True)
# ~ abserr = 0.0035 # absolute error of traveltimes (s)
# ~ srtData.set('err', np.ones(srtData.size()) * abserr)
# ~ srtData.set('t', srtData['t']-.02)
# ~ srtData.markInvalid(srtData['t'] <= 0.)
# ~ srtData.removeInvalid()

# ~ # create TravelTimeManger object
# ~ mgr = TravelTimeManager(srtData)

# ~ # create mesh
# ~ plc = mt.createParaMeshPLC(srtData.sensorPositions(),
                           # ~ paraDX=0.25, # number of nodes between stations (sensor positions)
                           # ~ paraDepth=15, # model depth (m) (0 --> estimated)
                           # ~ paraBoundary=2, # margin for parameter domain 
                                           # ~ # (sensor spacing)
                           # ~ paraMaxCellSize=.5, # maximum cell area (m^2)
                           # ~ boundary=0) # boundary around the parameter domain
                                       # ~ # no need to worry about this here
# ~ mesh = mt.createMesh(plc,
                     # ~ quality=32, # 2D triangle quality
                                 # ~ # recommended values: 30 - 34
                     # ~ smooth=(1, 10)) # smoothing alogrithm - no need to change this

# ~ # set the mesh
# ~ mgr.setMesh(mesh,
            # ~ secNodes=3) # secondary nodes to refine the mesh to
                        # ~ # permit a more accurate ray tracing

# ~ # define the inversion parameters
# ~ lam = 10 # use different values
# ~ zWeight = 1 # isotropic smoothing
                        
# ~ # run inversion
# ~ vest = mgr.invert(vTop=500, # velocity at the top of the initial model (m/s)
                  # ~ vBottom=5000, # velocity at the bottom of the initial model (m/s)
                  # ~ lam=lam, # regularization parameter
                  # ~ zWeight=zWeight, # 1 - same smoothing in x and z direction
                                   # ~ # 0.25 - four times stronger smoothing in 
                                   # ~ #        in horizontal direction
                  # ~ robustData=False, verbose=True) # use iterative re-weighting

# ~ # Consider the ray coverage in the plot, i.e.,
# ~ # blank areas in the model that should not be
# ~ # interpreted.
# ~ scov = mgr.standardizedCoverage()
# ~ vest = np.ma.masked_where(scov == 0, vest)

# ~ # plot the inversion result
# ~ cMin = 750
# ~ cMax = 4000

# ~ fig, ax = pg.plt.subplots(1, figsize=(8, 5))

# ~ _, cb = pg.show(mesh, # mesh
                # ~ vest, # resolved velocities
                # ~ coverage=scov, # blank areas with low coverage
                # cMap='inferno_r', # colormap
                # ~ cMap='Spectral', # colormap
                # cMap='plasma', # colormap
                # ~ logScale=False, # illustrate log-scaled velocities
                # ~ cMin=cMin, # lower limit of colormap
                # ~ cMax=cMax, # upper limit of colormap
                # ~ colorBar=True, # show colorbar
                # ~ orientation='horizontal', # orientation of the colorbar
                # ~ ax=ax, # axis object to plot in
                # ~ fitView=True, # automatic aspect ratio
                # ~ showBoundary=False) # plot mesh boundaries

# ~ pg.show(geom, ax=ax, fillRegion=False, regionMarker=False)

# ~ # format plot
# ~ ax.set_xlabel('x (m)')
# ~ ax.set_ylabel('z (m)')
# ~ ax.spines['top'].set_visible(False)
# ~ ax.spines['right'].set_visible(False)
# ~ cb.set_label('v$_p$ ($\mathregular{m\/s^{-1}}$)', fontsize=10)

# ~ pg.plt.show()
