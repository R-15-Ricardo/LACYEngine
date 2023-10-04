from flask import Flask
from flask import request as flaskreq
from flask_restful import Api, Resource, reqparse
import requests
from flask_cors import CORS
from fastkml import KML
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon
import pandas as pd
import numpy as np

from pykrige.ok import OrdinaryKriging
from pyproj import Transformer

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import networkx as nx

from Modules.HUH import helpers

#hola

app = Flask(__name__)
cors = CORS(app)
api = Api(app)

def initialize_engine():
    # ----------- LOAD SAMPLE POINT VARIABLE DATA -----------

    # Read the shapefile
    river_data = gpd.read_file('data/santiago/rio_santiago_prueba.shp')
    river_data = river_data.to_crs("epsg:4326")

    # Read the CSV file
    measured_variables = pd.read_csv('data/rizo20092013update.csv')
    measured_variables = measured_variables.drop('Unnamed: 0', axis=1)

    # Read the KML file
    sampling_stations = helpers.kml_to_geodataframe('data/muestreo.csv.kml')
    sampling_stations = sampling_stations[:10].copy()
    sampling_stations['name'] = [f'RS-{i:02d}' for i in range(1,11)]

    # ----------- LOAD SAMPLE POINT GEOGRAPHIC DATA -----------

    # Extract the coordinates of the sampling stations
    coordinates = sampling_stations.geometry.apply(lambda point: [point.y, point.x]).tolist()

    # Get chosen variable df
    variable_name = 'Coliformes totales'
    var = measured_variables[['fecha', 'idPuntoMuestreo', variable_name]]

    # Apply to sampling stations
    values = []
    for i in range(1,11):
        val = (var[var['idPuntoMuestreo'] == np.float64(i)].dropna())['Coliformes totales'].astype(np.float64).values.mean()
        values.append(val)

    # ----------- GENERATE FULL AREA MESH GRID -----------

    # Prepare the data for Kriging
    lats, lons = zip(*coordinates)

    # Convert lat, lon to UTM
    lats_utm, lons_utm = helpers.latlon_to_utm(lats, lons)

    # Calculate the number of grid cells in UTM coordinates
    cell_size = 500  # 5000 meters
    lon_range_utm = max(lons_utm) - min(lons_utm)
    lat_range_utm = max(lats_utm) - min(lats_utm)

    lon_cells = int(np.ceil(lon_range_utm / cell_size))
    lat_cells = int(np.ceil(lat_range_utm / cell_size))

    grid_lon_utm_ls = np.linspace(min(lons_utm), min(lons_utm) + lon_cells * cell_size, num=lon_cells)
    grid_lat_utm_ls = np.linspace(min(lats_utm), min(lats_utm) + lat_cells * cell_size, num=lat_cells)

    grid_lon_utm, grid_lat_utm = np.meshgrid(grid_lon_utm_ls, grid_lat_utm_ls)
    # Convert the UTM grid back to lat, lon
    grid_lat, grid_lon = helpers.utm_to_latlon(grid_lat_utm, grid_lon_utm)

    # Grid dimension
    geo_grid_shape = grid_lat.T.shape

    # ----------- GENERATE MESH GRID MASK (BOUND CONDITIONS) -----------

    # Obtain river geodata (lat,lon)
    river_lats = []
    river_lons = []

    for line in list(river_data.geometry[0].geoms)[::-1]:
        for point in line.coords:
            river_lats.append(point[1])
            river_lons.append(point[0])

    start = 14000
    stop = -150

    river_lats = river_lats[start:stop]
    river_lons = river_lons[start:stop]

    # Generate river mask
    # Create a LineString from river_lons and river_lats
    river_line = LineString(list(zip(river_lons, river_lats)))

    # Create a mask for the meshgrid based on river_line
    mask = np.zeros(geo_grid_shape, dtype=bool)

    for i in range(geo_grid_shape[0]-1):
        for j in range(geo_grid_shape[1]-1):
            cell_polygon = Polygon([
                (grid_lon.T[i, j], grid_lat.T[i, j]),
                (grid_lon.T[i + 1, j], grid_lat.T[i + 1, j]),
                (grid_lon.T[i + 1, j + 1], grid_lat.T[i + 1, j + 1]),
                (grid_lon.T[i, j + 1], grid_lat.T[i, j + 1])
            ])

            mask[i, j] = not river_line.intersects(cell_polygon)

    # Handle last row
    for j in range(geo_grid_shape[1] - 1):
        cell_polygon = Polygon([
            (grid_lon.T[-2, j], grid_lat.T[-2, j]),
            (grid_lon.T[-1, j], grid_lat.T[-1, j]),
            (grid_lon.T[-1, j + 1], grid_lat.T[-1, j + 1]),
            (grid_lon.T[-2, j + 1], grid_lat.T[-2, j + 1])
        ])
        mask[-1, j] = not river_line.intersects(cell_polygon)

    # Handle last column
    for i in range(geo_grid_shape[0] - 1):
        cell_polygon = Polygon([
            (grid_lon.T[i, -2], grid_lat.T[i, -2]),
            (grid_lon.T[i + 1, -2], grid_lat.T[i + 1, -2]),
            (grid_lon.T[i + 1, -1], grid_lat.T[i + 1, -1]),
            (grid_lon.T[i, -1], grid_lat.T[i, -1])
        ])
        mask[i, -1] = not river_line.intersects(cell_polygon)

    # Handle bottom right corner cell
    cell_polygon = Polygon([
        (grid_lon.T[-2, -2], grid_lat.T[-2, -2]),
        (grid_lon.T[-1, -2], grid_lat.T[-1, -2]),
        (grid_lon.T[-1, -1], grid_lat.T[-1, -1]),
        (grid_lon.T[-2, -1], grid_lat.T[-2, -1])
    ])
    mask[-1, -1] = not river_line.intersects(cell_polygon)

    # ----------- GENERATE FLOW GRAPH -----------

    # Create an empty graph
    G = nx.DiGraph()

    # Add nodes to the graph based on the mask
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if not mask[i, j]:
                node = (i, j)
                G.add_node(node)

    # Find the southernmost grid cell
    south_node = None
    min_lat = np.inf

    for node in G.nodes:
        i, j = node
        lat = grid_lat.T[i, j]
        if lat < min_lat:
            min_lat = lat
            south_node = node

    # Perform DFS traversal from the southernmost grid cell
    visited = set()
    helpers.dfs(south_node, visited, G)

    monitoring_station = []
    for lat, lon in zip(lats, lons):
        closest_node = None
        min_distance = np.inf
        
        for node in G.nodes:
            i, j = node
            node_lat, node_lon = grid_lat.T[i, j], grid_lon.T[i, j]
            distance = np.sqrt((node_lat - lat)**2 + (node_lon - lon)**2)
            
            if distance < min_distance:
                min_distance = distance
                closest_node = node
        
        monitoring_station.append(closest_node)

    pos = {node: (grid_lon.T[node], grid_lat.T[node]) for node in G.nodes}  # Define positions based on lon and lat

    return values, lons_utm, lats_utm, grid_lon_utm_ls, grid_lat_utm_ls, mask, grid_lat, grid_lon, monitoring_station, G
    

class LACYPort(Resource):
    def get(self):
        values, lons_utm, lats_utm, grid_lon_utm_ls, grid_lat_utm_ls, mask, grid_lat, grid_lon, monitoring_station, G = app.config['inference_settings']
        n_steps = int(flaskreq.args.get('steps'))
        now_values = values

        time_k = []

        for _ in range(n_steps):
            # -------- Generate OK model for initial values --------
            print('Calculating OK')
            ok = OrdinaryKriging(lons_utm, lats_utm, now_values, variogram_model='linear')
            # Estimate values and variance on the grid in UTM coordinates
            z_estimated, _ = ok.execute('grid', grid_lon_utm_ls, grid_lat_utm_ls)
            z_masked = np.ma.array(z_estimated.T, mask=mask)
            val = {node: z_masked[node] for node in G.nodes}

            # -------- Plot values for current state --------

            density_lat = []
            density_lon = []
            density_z = []

            for i in range(z_masked.shape[0]):
                for j in range(z_masked.shape[1]):
                    if not mask[i, j]:
                        density_lat.append(grid_lat.T[i,j])
                        density_lon.append(grid_lon.T[i,j])
                        density_z.append(z_masked.data[i,j])

            time_k.append((density_lat,density_lon,density_z))

            #-------- Generate new state --------
            
            print('generating timestep')
            for _ in range(50):
                val = helpers.next_state(G,val)

            now_values = [val[node]*0.75 for node in monitoring_station]
        
        return time_k

api.add_resource(LACYPort, "/")

# with app.app_context():
# @app.before_first_request
# def preload_variables():
    # app.config['inference_settings'] = initialize_engine()

if __name__ == "__main__":
    with app.app_context():
        app.config['inference_settings'] = initialize_engine()
        
    app.run(port=5000, debug=False)
    