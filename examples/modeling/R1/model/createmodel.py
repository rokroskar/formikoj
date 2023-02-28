import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg
import pygimli.meshtools as mt

from formikoj import SeismicWaveformModeler

# create the model geometry

top = mt.createPolygon([[0, 0], [94, 0], 
                        [94, -3.5], [72, -3.5], 
                        [20, -2], [0, -2]],
                       isClosed=True, marker=1, area=0.1)

infill = mt.createPolygon([[20, -2], [72, -3.5], 
                           [70, -6], [22, -6]],
                          isClosed=True, marker=2, area=0.1)

bottom = mt.createPolygon([[0, -2], [20, -2], 
                           [22, -6], [70, -6], 
                           [72, -3.5], [94, -3.5], 
                           [94, -10], [0, -10]],
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
mesh.save('mesh.bms')

# define velocity values for the different model regions
velmap = [[1, 750],
          [2, 2500],
          [3, 4000]]
          
# show the final parameterized model
pg.show(mesh, data=velmap, label=pg.unit('vel'), showMesh=True)
plt.show()
