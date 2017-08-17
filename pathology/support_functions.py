from shapely.geometry import Polygon, Point
import os, re, math
import numpy as np
from scipy.misc import imread


# takes in coordinates of a tile, its size, and a scaled polygon,
# returns the ratio of malignant to beging points in this tile
def malingnancy(tile, tile_size, polygon):
    tile_x, tile_y = tile[0]*tile_size, tile[1]*tile_size
    malignant, benign = 0, 0
    poly = Polygon(polygon)
    for x in range(tile_x, tile_x + tile_size):
        for y in range(tile_y, tile_y + tile_size):
            p = Point(x,y)
            if (poly.contains(p)):
                malignant+=1
            else:
                benign+=1
    if benign == 0: #patch is entirely malignant
        return 65536.0
    return malignant/benign

# takes in coordinates of a tile, its size, and list of all scaled regions,
# (regions as tuples: (bounds, vertex coordinates));
# returns the ratio of all malignant to bening points in this tile
# the difference is that this checks for all possible malignancies, not just
# for the current region, like previous function
def total_malignancy(tile, tile_size, scaled_regions):
    tile_x, tile_y = tile[0]*tile_size, tile[1]*tile_size
    malignant, benign = 0, 0
    fully_benign = True
    for region_bounds, region_coordinates in scaled_regions:
        if tile_in_region(tile, tile_size, region_bounds):
            fully_benign = False
            poly = Polygon(region_coordinates)
            for x in range(tile_x, tile_x + tile_size):
                for y in range(tile_y, tile_y + tile_size):
                    p = Point(x,y)
                    if (poly.contains(p)):
                        malignant+=1
                    else:
                        benign+=1
    if benign == 0: #patch is entirely malignant
        return 65536.0
    if fully_benign:
        return 0.0
    return malignant/benign

# returns the ratio of white to non-white points here
def whitespace(tile_size, tile):
    count = 0
    picture = np.ndarray.flatten(tile)
    for pixel in picture:
        if pixel > 230:
            count += 1
    return count / len(picture)

#finds a 'bounding box' for a given region
def region_bounds(coordinates):
    x_coord, y_coord = [], []
    for x,y in coordinates:
        x_coord.append(x)
        y_coord.append(y)
    return min(x_coord), max(x_coord), min(y_coord), max(y_coord)

#goes to given directory or creates it if it's not there
def make_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

#parses the dzi string for image dimensions
def get_image_dimensions(dzi):
    p = re.compile('Height=\"(.*?)\"')
    image_height = int(p.search(dzi).group(1))
    p = re.compile('Width=\"(.*?)\"')
    image_width = int(p.search(dzi).group(1))
    return image_height, image_width

def scale_region(highest_zoom, zoom_level, region_bounds, region):
    scaling_factor = math.pow(2,(highest_zoom-zoom_level-1))
    scaled_region_bounds = [n / scaling_factor for n in region_bounds]
    scaled_region = [(n / scaling_factor, m / scaling_factor)
    for (n,m) in region]
    return scaled_region_bounds, scaled_region

def tile_in_region(tile, tile_size, scaled_region_bounds):
    tile_x, tile_y = tile[0]*tile_size, tile[1]*tile_size
    min_x, max_x, min_y, max_y = scaled_region_bounds
    x_in_region = ((tile_x > min_x and tile_x + tile_size < max_x) or
    (tile_x < min_x and tile_x + tile_size > min_x) or
    (tile_x < max_x and tile_x + tile_size > max_x) or
    (tile_x < min_x and tile_x + tile_size > max_x))
    y_in_region = ((tile_y > min_y and tile_y + tile_size < max_y) or
    (tile_y < min_y and tile_y + tile_size > min_y) or
    (tile_y < max_y and tile_y + tile_size > max_y) or
    (tile_y < min_y and tile_y + tile_size > max_y))
    # if (x_in_region and y_in_region):
    #     print("Tile: []%d, %d, %d, %d] IS in region: " %
    #     (tile_x, tile_x+tile_size, tile_y, tile_y+tile_size), scaled_region_bounds)
    # else:
    #     print("Tile: []%d, %d, %d, %d] IS NOT in region: " %
    #     (tile_x, tile_x+tile_size, tile_y, tile_y+tile_size), scaled_region_bounds)
    return x_in_region and y_in_region
