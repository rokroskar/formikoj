# -*- coding: utf-8 -*-

# Import Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import pygimli as pg
import pyvista as pv
from pyvista import examples
import pandas as pd
from scipy.interpolate import griddata

from datetime import datetime
# ~ import vtk

from plot_3D_tools import read_raster

# read the pck file
srtData = pg.DataContainer("../picks/picks_20220714T104939.pck", # change the name of your pick file
                           sensorTokens="s g",
                           removeInvalid=True)

sensor_coords = np.array(srtData.sensorPositions())
sgt = np.vstack((np.array(srtData['s']),
                 np.array(srtData['g']), 
                 np.array(srtData['t']))).T

# ~ vp = pv.read('3dsb_sn3.vtk')
# ~ slice001 = vp.slice_orthogonal()
# ~ slice002 = vp.slice(normal=[1, 1, 0])
# ~ slice003 = vp.slice(normal=[-1, 1, 0])

#%% Show 3D surface with threshold PLOT
# Import DEM
terrain = read_raster('dem_sonnblick_clipped_gkm31.tif')
# ~ glacier = read_raster('dem_sonnblick_clipped_glacier.tif')

# Get terrain
terrain = terrain.warp_by_scalar(factor=0)
# ~ glacier = glacier.warp_by_scalar(factor=0)

# Shift hut and terrain
# ~ glacier.translate([0,0,0.1])

# Get contour lines
lines = np.arange(3020, 3110, 10)
contours = terrain.contour(lines)

# Create a pyvista plotter
# ~ p = pv.Plotter(notebook=False, off_screen=True, window_size=(1900, 900),
p = pv.Plotter(notebook=False, off_screen=True, window_size=(1200, 800),
               )

# Add terrain
# ~ p.add_mesh(glacier,
            # ~ color='lightskyblue', show_edges=False,
            # ~ opacity=1)
p.add_mesh(terrain,
            color='grey', show_edges=False, edge_color='grey',
            opacity=.3, lighting=False)
# ~ p.add_mesh(vp)
sbar_args = {"n_labels": 2, "title": "Seismic Velocity (m/s)", "color": "k", "fmt": "%.0f"}
p.add_mesh(slice001, scalars="Velocity", clim=[500, 5000], scalar_bar_args=sbar_args, cmap="coolwarm_r") #,
           # ~ use_transparency=False, opacity="S_Coverage")
p.add_mesh(slice002, scalars="Velocity", clim=[500, 5000], show_scalar_bar=False, cmap="coolwarm_r") #,
           # ~ use_transparency=False, opacity="S_Coverage")
p.add_mesh(slice003, scalars="Velocity", clim=[500, 5000], show_scalar_bar=False, cmap="coolwarm_r") #,
           # ~ use_transparency=False, opacity="S_Coverage")
# ~ # Import electrodes
# ~ elec_nw = np.loadtxt('sbk_nordwand_borehole_shift_utm_newz.csv', skiprows=1, 
                     # ~ delimiter=',', usecols=(1,2,4))
# ~ elec_mo = np.loadtxt('sbk_gba_monitoring_utm.csv', skiprows=1, 
                     # ~ delimiter=',', usecols=(0,1,5))
# ~ elec_gl = np.loadtxt('sbk_line3_sip_utm.csv', skiprows=1, 
                     # ~ delimiter=',', usecols=(1,2,3))

# ~ # Shift electrodes
# ~ elec_nw[:,2] = elec_nw[:,2]+1
# ~ elec_nw[:32,1] = elec_nw[:32,1]+4
# ~ elec_mo[:,2] = elec_mo[:,2]+1
# ~ elec_gl[:,2] = elec_gl[:,2]+1


# ~ # Create polydata
# ~ elec_nw_pts = pv.PolyData(elec_nw)
# ~ elec_mo_pts = pv.PolyData(elec_mo)
# ~ elec_gl_pts = pv.PolyData(elec_gl)

# ~ # Plot electrodes
# ~ p.add_mesh(elec_nw_pts, color='black', point_size=5,
            # ~ render_points_as_spheres=False)
# ~ p.add_mesh(elec_mo_pts, color='black', point_size=5,
            # ~ render_points_as_spheres=False)
# ~ p.add_mesh(elec_gl_pts, color='black', point_size=5,
            # ~ render_points_as_spheres=False)

# # Add contour lines
# ~ p.add_mesh(contours, line_width=3, color='white')

# ~ # Add text
# ~ p.add_point_labels([elec_nw[20,0]+10, elec_nw[20,1], elec_nw[20,2]],
                    # ~ labels=['North wall'], show_points=False,
                    # ~ font_size=21)

# ~ p.add_point_labels([elec_nw[0,0], elec_nw[0,1], elec_nw[0,2]],
                    # ~ labels=['P4'], show_points=False,
                    # ~ font_size=16)

# ~ p.add_point_labels([elec_nw[63,0], elec_nw[63,1], elec_nw[63,2]],
                    # ~ labels=['P2'], show_points=False,
                    # ~ font_size=16)

# ~ p.add_point_labels([elec_gl[31,0], elec_gl[31,1], elec_gl[31,2]],
                    # ~ labels=['P3'], show_points=False,
                    # ~ font_size=16)

# ~ p.add_point_labels([elec_mo[0,0], elec_mo[0,1], elec_mo[0,2]],
                    # ~ labels=['P1'], show_points=False,
                    # ~ font_size=16)

# ~ p.add_point_labels([elec_nw[50,0], elec_nw[50,1], elec_nw[50,2]],
                    # ~ labels=['South slope'], show_points=False,
                    # ~ font_size=21)

# ~ p.add_point_labels([elec_gl[31,0]-30, elec_gl[31,1]-50, elec_gl[31,2]-10],
                    # ~ labels=['Glacier'], show_points=False,
                    # ~ font_size=21)

# ~ p.add_point_labels([elec_nw[32,0]+15, elec_nw[32,1]-15, elec_nw[32,2]+15],
                    # ~ labels=['Hoher Sonnblick (3106 m)'], show_points=False,
                    # ~ font_size=21)

# Camera
p.camera.azimuth = 225
p.camera.elevation = -15
p.camera.zoom(2.0)

# Show axes
# p.show_bounds(show_yaxis=True,
#               xlabel='x (m)',
#               ylabel='y (m)',
#               zlabel='z (m)',
#               color='black',
#               location='origin',
#               grid='back')

# Light
# light = pv.Light(intensity=0.7)
# light.set_direction_angle(30, 290)
# p.add_light(light)

# Set white background
p.set_background(color='white')

# Show plot
# p.show(full_screen=True,
#        #cpos=[-1, -0.1, 0.6],
#        screenshot = 'map_v01.png'
#        )

# Video
orbit_path = p.generate_orbital_path(n_points=100, shift=terrain.length/3)
p.open_gif('map_v01.gif')
p.orbit_on_path(orbit_path, write_frames=True, step=1)
p.close()


# plotter = pv.Plotter()
# plotter.add_mesh(elec_gl_pts, color='maroon', style='points',
#                  point_size=10.0)
# # plotter.add_point_labels([point_cloud.center,], ['Center',],
# #                          point_color='yellow', point_size=20)
# #plotter.show_grid()
# plotter.show()
