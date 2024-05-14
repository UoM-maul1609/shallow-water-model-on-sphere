"""
Module Name: create_map.py

Description:
Using the mpl_toolkits Basemap tool to create a 
Basemap object, draw meridians and parallels,
and return this object along with the x and y 
map coordinate matrices.
"""

import numpy as np
from mpl_toolkits.basemap import Basemap

def create_map_func(lons, lats, satellite_height=3000000, satellite_angle=90):
    """
    Draw an orthographic map projection with perspective of satellite 
    looking down at 50N, 100W.

    Parameters:
        lons (numpy.ndarray): 1D array of longitudes.
        lats (numpy.ndarray): 1D array of latitudes.
        satellite_height (float): The height we want the map to be
                                  displayed from, from Saturn's core.
        satellite_angle (float): Angle from which to view the Saturn
                                 map.

    Returns:
        tuple: A tuple containing the meshgrid x and y coordinates.
                x and y are matrices of the same size as data, containing
                the positions of the elements in the map coordinates.
               And also the Basemap object.
    """
    basemap = Basemap(
        satellite_height=satellite_height,
        projection="nsper",
        lat_0=satellite_angle,
        lon_0=-100,
        resolution="l",
    )

    # Draw the edge of the map projection region (the projection limb)
    # Draw lat/lon grid lines every 30 degrees.
    basemap.drawmeridians(np.arange(0, 360, 30))
    basemap.drawparallels(np.arange(-90, 90, 30))

    # Make sure that map boundary is not clipped by the axes
    circle = basemap.drawmapboundary(linewidth=0.75, color='k')
    circle.set_clip_on(False)

    # Create meshgrid
    (lo, la) = np.meshgrid(lons, lats)
    x, y = basemap(lo * 180.0 / np.pi, la * 180.0 / np.pi)

    return x, y, basemap
