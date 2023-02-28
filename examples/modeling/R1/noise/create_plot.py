import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

from obspy import read

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import TravelTimeManager

nstats = 48
spc = 2

stations = np.vstack((np.linspace(0, (nstats-1)*spc, nstats), 
                      np.zeros(nstats))).T

fig, ax = plt.subplots(4, 2, figsize=(8, 10), #sharex=True,
                       gridspec_kw={'width_ratios': [1., .02]})

# PLOT MODEL
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

cMin = 750
cMax = 4000

mesh = pg.load('./in/mesh.bms')
velmap = [[1, 750],
          [2, 2500],
          [3, 4000]]
          
pg.show(mesh, data=velmap,
        cMap='Spectral', # colormap
        logScale=False, # illustrate log-scaled velocities
        cMin=cMin, # lower limit of colormap
        cMax=cMax, # upper limit of colormap
        colorBar=False, # show colorbar
        orientation='vertical', # orientation of the colorbar
        fitView=False, # automatic aspect ratio
        showBoundary=False,
        ax=ax[0, 0])
        
cmap = mpl.cm.Spectral
norm = mpl.colors.Normalize(vmin=cMin, vmax=cMax)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax[0, 1], orientation='vertical', label='P-wave velocity (m/s)')
        
ax[0, 0].plot(stations[:23, 0], stations[:23, 1], ls='',
           marker='v', color='w', markersize=6, markeredgecolor='k')
ax[0, 0].plot(stations[23, 0], stations[23, 1],
           marker='v', color='w', markersize=6, markeredgecolor='k')
ax[0, 0].plot(stations[24:, 0], stations[24:, 1], ls='',
           marker='v', color='w', markersize=6, markeredgecolor='k')

ax[0, 0].set_yticklabels(np.arange(0, 20, 5)[::-1])

# PLOT WAVEFORM DATA
st = read('/media/matthias/Data/git/formikoj/examples/modeling/R1/noise/out/syn_data/data/Shot_1024.syn')

t = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate

# read the pck file
srtData = pg.DataContainer("./out/syn_data/syn_data_tt.pck", 
                           sensorTokens="s g",
                           removeInvalid=True)
                           
sgt = np.vstack((srtData['s'].copy(), srtData['g'].copy(), srtData['t'].copy())).T

sgt[:, 0] *= spc
sgt[:, 1] *= spc

sgt[sgt[:, 1]==4., 2] = np.nan
sgt[sgt[:, 1]==88., 2] = np.nan
sgt[sgt[:, 1]==70., 2] += .025

# plot seismograms
X, Y = [], []
fill = []
vd = np.empty((len(t), len(st)))
xtls = np.empty(len(st), dtype=int)
for i, tr in enumerate(st):
    # prepare trace data for plotting
    trd = tr.data.copy()
                
    trd = trd * 12 + i*spc
    
    # clip amplitudes
    trd[trd > .45 + i*spc] = .45 + i*spc; trd[trd < i*spc-.45] = i*spc-.45
    
    X.append(list(trd))
    Y.append(list(t))

    # populate variable density array
    vd[:, i] = trd - i*spc

# add seismograms to plot
seismograms = LineCollection((list(zip(x, y)) for x, y in zip(X, Y)),
                             colors='k', linestyles='-',
                             linewidths=.5)
ax[1, 0].add_collection(seismograms)

ax[1, 0].imshow(vd,
          extent=(-.5, i*spc + .5, t[-1], t[0]),
          cmap='seismic_r',
          aspect='auto',
          vmin=-1, vmax=1,
          interpolation='bicubic')

idx = sgt[:, 0] == 23*spc
ax[1, 0].scatter(sgt[idx, 1], sgt[idx, 2], s=30, marker='X', c='g',
                 zorder=100, edgecolor='k')

ax[1, 0].plot(stations[:23, 0], stations[:23, 1]-.005, ls='',
           marker='v', color='w', markersize=6, markeredgecolor='k', zorder=100)
ax[1, 0].plot(stations[23, 0], stations[23, 1]-.005,
           marker='*', color='w', markersize=8, markeredgecolor='k', zorder=100)
ax[1, 0].plot(stations[24:, 0], stations[24:, 1]-.005, ls='',
           marker='v', color='w', markersize=6, markeredgecolor='k', zorder=100)

ax[1, 0].set_ylim((.11, -.01))

# PLOT PSEUDOSECTION

cmap = mpl.cm.Spectral
norm = mpl.colors.Normalize(vmin=300, vmax=1500)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax[2, 1], orientation='vertical', label='Apparent velocity (m/s)')

mp = (sgt[:, 0] + sgt[:, 1])/2
d = np.abs(sgt[:, 1] - sgt[:, 0])
pd = np.round(d/5, 0)
av = np.round(d/sgt[:, 2], 1)

mult = np.vstack((mp, pd)).T

# find unique mp-pd pairs
u, idx, inv, cnt = np.unique(mult,
                             return_counts=True,
                             return_index=True,
                             return_inverse=True,
                             axis=0)

# plot pseudosection
for i in np.arange(np.max(cnt)) + 1:
    scat = ax[2, 0].scatter(mp[idx[cnt == i]],
                        pd[idx[cnt == i]],
                        c=av[idx[cnt==i]], s=40,
                        cmap=cmap,
                        norm=norm)
    if i > 1:
        # get idices of multiple pairs
        mi = [np.where(np.logical_and(
            mult[:,0]==mult[j, 0], \
            mult[:,1]==mult[j, 1]) == 1)[0][-1]\
                for j in idx[cnt == i]]
                
        ax[2, 0].scatter(mp[mi], 
                         pd[mi],
                     c=av[mi], s=40 / i ** 3,
                     cmap=cmap, 
                     norm=norm)

ax[2, 0].plot(stations[:, 0], stations[:, 1]-.75, ls='',
           marker='v', color='w', markersize=6, markeredgecolor='k')
ax[2, 0].invert_yaxis()

# PLOT INVERSION RESULT
# read the pck file
srtData = pg.DataContainer("./out/syn_data/syn_data_tt.pck", 
                           sensorTokens="s g",
                           removeInvalid=True)
abserr = 0.0003 # absolute error of traveltimes (s)
srtData.set('err', np.ones(srtData.size()) * abserr)
srtData.set('t', srtData['t']-.02)
srtData.markInvalid(srtData['t'] <= 0.)
srtData.removeInvalid()

# create TravelTimeManger object
mgr = TravelTimeManager(srtData)

# create mesh
plc = mt.createParaMeshPLC(srtData.sensorPositions(),
                           paraDX=0.25, # number of nodes between stations (sensor positions)
                           paraDepth=15, # model depth (m) (0 --> estimated)
                           paraBoundary=2, # margin for parameter domain 
                                           # (sensor spacing)
                           paraMaxCellSize=.25, # maximum cell area (m^2)
                           boundary=0) # boundary around the parameter domain
                                       # no need to worry about this here
invmesh = mt.createMesh(plc,
                     quality=32, # 2D triangle quality
                                 # recommended values: 30 - 34
                     smooth=(1, 10)) # smoothing alogrithm - no need to change this

# set the mesh
mgr.setMesh(invmesh,
            secNodes=3) # secondary nodes to refine the mesh to
                        # permit a more accurate ray tracing

# define the inversion parameters
lam = 10 # use different values
zWeight = .5 # isotropic smoothing
                        
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
_, cb = pg.show(invmesh, # mesh
                vest, # resolved velocities
                coverage=scov, # blank areas with low coverage
                cMap='Spectral', # colormap
                logScale=False, # illustrate log-scaled velocities
                cMin=cMin, # lower limit of colormap
                cMax=cMax, # upper limit of colormap
                colorBar=False, # show colorbar
                orientation='vertical', # orientation of the colorbar
                ax=ax[3, 0], # axis object to plot in
                fitView=False, # automatic aspect ratio
                showBoundary=False) # plot mesh boundaries

pg.show(geom, ax=ax[3, 0], fillRegion=False, regionMarker=False, fitView=False)

ax[3, 0].plot(stations[:23, 0], stations[:23, 1], ls='',
           marker='v', color='w', markersize=6, markeredgecolor='k')
ax[3, 0].plot(stations[23, 0], stations[23, 1],
           marker='v', color='w', markersize=6, markeredgecolor='k')
ax[3, 0].plot(stations[24:, 0], stations[24:, 1], ls='',
           marker='v', color='w', markersize=6, markeredgecolor='k')

ax[0, 0].set_yticklabels(np.arange(0, 20, 5)[::-1])

cmap = mpl.cm.Spectral
norm = mpl.colors.Normalize(vmin=cMin, vmax=cMax)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax[3, 1], orientation='vertical', label='P-wave velocity (m/s)')


ax[0, 0].set_ylim((-15, 1))
ax[3, 0].set_ylim((-15, 1))

for a in ax[:, 0]:
    a.set_xlim((-2.5, 94+2.5))
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
    a.set_xlabel('Distance (m)')

ax[0, 0].set_ylabel('Depth (m)')
ax[1, 0].set_ylabel('Time (s)')
ax[2, 0].set_ylabel('Pseudodepth (m)')
ax[3, 0].set_ylabel('Depth (m)')

ax[1, 1].set_visible(False)

ax[0, 0].set_title('a)', loc='left')
ax[1, 0].set_title('b)', loc='left')
ax[2, 0].set_title('c)', loc='left')
ax[3, 0].set_title('d)', loc='left')

plt.tight_layout()
# ~ plt.show()
# ~ plt.savefig('synthetic_data_example.png', dpi=300)
plt.savefig('synthetic_data_example.pdf')
