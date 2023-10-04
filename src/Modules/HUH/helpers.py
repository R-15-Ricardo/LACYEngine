from fastkml import KML
import geopandas as gpd
from shapely.geometry import Point
from pyproj import Transformer

def kml_to_geodataframe(kml_file):
    with open(kml_file, 'rb') as f:
        kml_data = f.read()

    k = KML()
    k.from_string(kml_data)

    features = list(k.features())[0].features()

    # Extract the coordinates and properties from the KML features
    coordinates = []
    properties = []

    for feature in features:
        coords = (feature.geometry.x, feature.geometry.y)
        props = {'name': feature.name}
        coordinates.append(coords)
        properties.append(props)

    # Create a GeoDataFrame from the extracted data
    geometry = [Point(coord) for coord in coordinates]
    gdf = gpd.GeoDataFrame(properties, geometry=geometry)

    return gdf

# Function to convert lat, lon to UTM
def latlon_to_utm(lats, lons):
    transformer = Transformer.from_crs("epsg:4326", "epsg:32613")  # Updated to EPSG:32613
    return transformer.transform(lats, lons)

# Function to convert UTM to lat, lon
def utm_to_latlon(xs, ys):
    transformer = Transformer.from_crs("epsg:32613", "epsg:4326")  # Updated to EPSG:32613 EPSG:26913
    return transformer.transform(xs, ys)

def dfs(node, visited, G):
    i, j = node
    visited.add(node)
    
    for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        neighbor_node = (i + di, j + dj)
        if neighbor_node in G and neighbor_node not in visited:
            G.add_edge(node, neighbor_node)
            dfs(neighbor_node, visited, G)

def next_state(G, val):
    new_val = {node: 0 for node in G.nodes}
    
    for node in G.nodes:
        out_edges = list(G.out_edges(node))
        if not out_edges:
            new_val[node] += val[node]  # Keep the value if there are no outgoing edges
        else:
            # Distribute the value equally among the adjacent nodes
            distributed_value = val[node] / len(out_edges)
            for _, neighbor_node in out_edges:
                new_val[neighbor_node] += distributed_value
                
    return new_val
