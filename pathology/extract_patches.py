import xml.etree.ElementTree as ET
from pathology_config import PATHOLOGYconfig
import openslide as ops
from openslide import ImageSlide, open_slide
from openslide.deepzoom import DeepZoomGenerator
from PIL import Image
from shapely import geometry
from matplotlib import pyplot as plt
from scipy.misc import imread
from scipy import ndimage
import numpy as np
import math
import support_functions as sf
from random import randint

#get all directory parameters from config
config = PATHOLOGYconfig()
input_xml = config.input_xml
home_folder = config.home_folder
patch_folder = home_folder + 'patches/'
sf.make_dir(patch_folder)

#get all patch parameters from config
patch_size = config.patch_size
patch_format = config.patch_format
colors = config.annotation_colors
number_benign_patches = config.number_benign_patches
combine_annotations = config.benign_combine_annotations
benign_maximum_malignancy = config.benign_maximum_malignancy

#print a report on what extract_patches will be doing?
verbose = config.verbosity

#parse the annotation .xml file
annotations = []
xml_tree = ET.parse(input_xml)
xml_root = xml_tree.getroot()
for annotation in xml_root.iter('Annotation'):
    region_coordinates, region_zoom_levels = [] , []
    annotation_id = int(annotation.get('Id'))
    for region in annotation.iter('Region'):
        region_id = int(region.get('Id'))
        region_zoom_level = float(region.get('Zoom'))
        coordinates = []
        for index, vertex in enumerate(region.iter('Vertex')):
            x, y = float(vertex.get("X")), float(vertex.get("Y"))
            coordinates.append((x,y))
        region_coordinates.append(coordinates)
        region_zoom_levels.append(region_zoom_level)

    annotations.append((region_coordinates,region_zoom_levels))
print("XML processed.")

#process svs
input_svs = config.svs
slide = ops.OpenSlide(input_svs)
zoom_gen = ops.deepzoom.DeepZoomGenerator(slide, tile_size=patch_size)
dzi = zoom_gen.get_dzi('jpeg')
image_height, image_width = sf.get_image_dimensions(dzi)

#create an output folder for the given slide
slide_name_split = input_svs.split('/')
short_slide_name = slide_name_split[len(slide_name_split)-1].split('.')[0]
slide_output_folder = patch_folder + short_slide_name + '/'
sf.make_dir(slide_output_folder)
output_malignant = slide_output_folder + "malignant/"
output_benign = slide_output_folder + "benign/"
output_random = slide_output_folder + "random/"
sf.make_dir(output_malignant)
sf.make_dir(output_benign)
sf.make_dir(output_random)

#create a thumbnail of the image with annotation
tile_size = max(image_height,image_width)/40
thumbnail_gen = ops.deepzoom.DeepZoomGenerator(slide, tile_size=1791.975)
thumbnail = thumbnail_gen.get_tile(11,(0,0))
img = np.flipud(thumbnail)
fig = plt.figure(1, figsize=(18,6), dpi=90)
a = plt.gca()
ax = fig.add_subplot(111)
ax.imshow(img, extent=[0,71679,0,23739])
ax.set_title('Malignant regions')
for annotation in range(len(annotations)):
    for region in range(len(annotations[annotation][0])):
        poly = geometry.Polygon(
        annotations[annotation][0][region])
        x, y = poly.exterior.xy
        ax.plot(x, y, color=colors[annotation])
        ax.axis([0,image_width,0,image_height])
a.invert_yaxis()
plt.savefig("Malignant_regions_plot", format='png')
print("Created annotated thumbnail at: %s." % slide_output_folder + "Malignant_regions_plot")
# plt.show()

#create patches of malignant tissue
#calculate reasonable zoom levels:
highest_zoom = zoom_gen.level_count #full resolution
#lowest_zoom: zoom level where whole image is less than one patch
lowest_zoom = math.floor(highest_zoom - (math.log(
min(image_width,image_height)/patch_size,2)))

#scale all regions for all zoom levels
regions_by_zoom_level_by_annotation = []
for ann_ind, annotation in enumerate(annotations):
    regions_by_zoom_level = []
    #scale all regions by zoom level
    for zoom_level in range(lowest_zoom, highest_zoom):
        regions_at_this_zoom_level = []
        for reg_ind, region in enumerate(annotations[ann_ind][0]):
            unscaled_bounds, unscaled_region = sf.region_bounds(
            annotations[ann_ind][0][reg_ind]), annotations[ann_ind][0][reg_ind]

            regions_at_this_zoom_level.append(sf.scale_region(
            highest_zoom, zoom_level, unscaled_bounds, unscaled_region))
        regions_by_zoom_level.append(regions_at_this_zoom_level)
    regions_by_zoom_level_by_annotation.append(regions_by_zoom_level)

#save malignant patches
def malignancy_ratio_check(verbose, zoom_generator, zoom_level, tile_count, patch_size, scaled_region_bounds, scaled_region, ann_ind, reg_ind, patch_format, folder, args):
    for x in range(tile_count[0]):
        for y in range(tile_count[1]):
            verbose_report = ""
            if sf.tile_in_region((x,y), patch_size, scaled_region_bounds):
                malignancy = sf.malingnancy((x,y), patch_size, scaled_region)
                if malignancy > 1:
                    patch = zoom_generator.get_tile(zoom_level, (x,y))
                    patch.save(folder+"ann_%d_reg_%d_zoom_%d_coord_%d_%d_malign_%d" % (
                    ann_ind, reg_ind, zoom_level, x, y, int(malignancy)), patch_format)
                    verbose_report = "--> saving this patch as malignant"
                if verbose:
                    print("Malignant/benign pts in the potentially " +
                    "malign. patch (%d, %d) : %f. %s" % (x,y, malignancy, verbose_report))

#save benign patches
def benign_ratio_check(verbose, zoom_generator, zoom_level, tile_count, patch_size, scaled_region_bounds, scaled_regions, ann_ind, reg_ind, patch_format, folder, args):
    verbose_report = ""
    number_benign_patches = args['number_benign_patches']
    all_scaled_regions, scaled_regions = args['scaled_regions'], []
    patch_whitespace_percent = args['benign_whitespace']
    count, iterations = 0, 0
    if combine_annotations:
        for ann_ind, annotation in enumerate(annotations):
            scaled_regions += all_scaled_regions[ann_ind][zoom_level-lowest_zoom]
    else:
        scaled_regions = all_scaled_regions[ann_ind][zoom_level-lowest_zoom]
    while count <= number_benign_patches and iterations <= tile_count[0]*tile_count[1]:
        iterations+=1
        a, b = randint(0,tile_count[0]-1), randint(0,tile_count[1]-1)
        total_malignancy = sf.total_malignancy((a,b), patch_size, scaled_regions)
        if total_malignancy < benign_maximum_malignancy:
            patch = zoom_generator.get_tile(zoom_level, (a,b))
            whitespace = sf.whitespace(patch_size, np.array(patch))
            if whitespace < patch_whitespace_percent:
                count+=1
                patch.save(folder+"ann_%d_reg_%d_zoom_%d_coord_%d_%d_malign_%d" % (
                ann_ind, reg_ind, zoom_level, a, b, int(total_malignancy)), patch_format)
                verbose_report = "--> saving this patch as benign"
            if verbose:
                print("Malignant/benign pts in the potentially " +
                "benign. patch (%d, %d) : %f, has %f whitespace ratio %s" % (a,b, total_malignancy, whitespace, verbose_report))
                verbose_report=""

#save benign patches
def random_check(verbose, zoom_generator, zoom_level, tile_count, patch_size, scaled_region_bounds, scaled_regions, ann_ind, reg_ind, patch_format, folder, args):
    number_benign_patches = args['number_random_patches']
    all_scaled_regions, scaled_regions = args['scaled_regions'], []
    patch_whitespace_percent = args['random_whitespace']
    count, iterations = 0, 0
    if combine_annotations:
        for ann_ind, annotation in enumerate(annotations):
            scaled_regions += all_scaled_regions[ann_ind][zoom_level-lowest_zoom]
    else:
        scaled_regions = all_scaled_regions[ann_ind][zoom_level-lowest_zoom]
    while count <= number_benign_patches and iterations <= tile_count[0]*tile_count[1]:
        iterations+=1
        a, b = randint(0,tile_count[0]-1), randint(0,tile_count[1]-1)
        patch = zoom_generator.get_tile(zoom_level, (a,b))
        whitespace = sf.whitespace(patch_size, np.array(patch))
        if whitespace < patch_whitespace_percent:
            count+=1
            total_malignancy = sf.total_malignancy((a,b), patch_size, scaled_regions)
            patch.save(folder+"ann_%d_reg_%d_zoom_%d_coord_%d_%d_malign_%d" % (
            ann_ind, reg_ind, zoom_level, a, b, int(total_malignancy)), patch_format)
    if verbose:
            print("Saved %d random patches (looking for less than %d percent whitespace)." % (count, patch_whitespace_percent*100))

def create_patches (malignancy_folder, checking_function, checking_func_args):

    for ann_ind, annotation in enumerate(annotations):
        annotation_folder = malignancy_folder + 'annotation_%d/' % ann_ind
        sf.make_dir(annotation_folder)

        for reg_ind, region in enumerate(annotations[ann_ind][0]):
            region_folder = annotation_folder + 'region_%d/' % reg_ind
            sf.make_dir(region_folder)
            region_bounds = sf.region_bounds(annotations[ann_ind][0][reg_ind])

            for zoom_level in range(lowest_zoom, highest_zoom):
                if verbose:
                    print("============= ZOOM_LEVEL: %d " % zoom_level +
                    " ============= REGION: %d" % reg_ind +
                    " ============= ANNOTATION: %d =============" % ann_ind)

                zoom_folder = region_folder+ 'zoom_level%d/' % zoom_level
                sf.make_dir(zoom_folder)

                # scale coordinates of the current region for the zoom level
                # if we're getting malignant patches, we need just the current region
                # if checking_function == malignancy_ratio_check:
                zoom_index = zoom_level-lowest_zoom
                scaled_region_bounds, scaled_region = regions_by_zoom_level_by_annotation[ann_ind][zoom_index][reg_ind]

                # find out what tiles are involved:
                tile_count = zoom_gen.level_tiles[zoom_level]

                # save pathes that satisfy the checking_function
                checking_function(verbose, zoom_gen, zoom_level, tile_count, patch_size, scaled_region_bounds, scaled_region,
                ann_ind, reg_ind, patch_format, zoom_folder, checking_func_args)

#create malignant patches
create_patches(output_malignant, malignancy_ratio_check, {})

#create benign patches
benign_ratio_check_arguments = {}
benign_ratio_check_arguments['benign_whitespace'] = config.benign_whitespace
benign_ratio_check_arguments['scaled_regions'] = regions_by_zoom_level_by_annotation
benign_ratio_check_arguments['number_benign_patches'] = number_benign_patches
create_patches(output_benign, benign_ratio_check, benign_ratio_check_arguments)

#create random patches
random_check_arguments = {}
random_check_arguments['number_random_patches'] = config.number_random_patches
random_check_arguments['random_whitespace'] = config.random_whitespace
random_check_arguments['scaled_regions'] = regions_by_zoom_level_by_annotation
create_patches(output_random, random_check, random_check_arguments)
