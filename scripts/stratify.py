from PIL import Image
import pickle
from collections import defaultdict
from multiprocessing import Pool
import os
import glob
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score
import numpy as np
import pandas as pd
import statistics
from scipy.spatial import distance
pd.options.mode.chained_assignment = None


cd45_color = '#F8766D'
ck_color = '#00BA38'
cd68_color = '#619CFF'
geomx_images = "/home/dimitris/Workspace/crc/data/images_RT016.22"
annotations_file = '/home/dimitris/Workspace/crc/data/geomx_annotations.csv'

def sort_key(item):
    try:
        return -item['scores']['CD68']['ck_score']
    except KeyError:
        return float('-inf')

def combine_images(columns, space, images):
    rows = len(images) // columns
    if len(images) % columns:
        rows +=  1
    width_max = max([Image.open(image).width for image in images])
    height_max = max([Image.open(image).height for image in images])
    background_width = width_max * columns + (space * (columns -  1))
    background_height = height_max * rows + (space * (rows -  1))
    background = Image.new('RGBA', (background_width, background_height), (255,  255,  255,  255))
    x =  0
    y =  0
    for i, image in enumerate(images):
        img = Image.open(image)
        x_offset = int((width_max - img.width) /  2)
        y_offset = int((height_max - img.height) /  2)
        background.paste(img, (x + x_offset, y + y_offset))
        x += width_max + space
        if (i +  1) % columns ==  0:
            y += height_max + space
            x =  0
    background.save(geomx_images  + '/combined_image.png')

def get_consecutive_integer_series(integer_list):
    integer_list = sorted(integer_list)
    start_item = integer_list[0]
    end_item = integer_list[-1]
    a = set(integer_list) # Set a
    b = range(start_item, end_item+1)
    # Pick items that are not in range.
    c = set(b) - a # Set operation b-a
    li = []
    start = 0
    for i in sorted(c):
        end = b.index(i) # Get end point of the list slicing
        li.append(b[start:end]) # Slice list using values
        start = end + 1 # Increment the start point for next slicing
    li.append(b[start:]) # Add the last series
    results = []
    for sliced_list in li:
        if not sliced_list or len(sliced_list) == 1:
            continue
        else:
            results.append((sliced_list[0], sliced_list[-1]))

    if len(results) > 0:
        return results
    else:
        return False

# Function to check if a group of y values could form a vertical line
def could_form_vertical_line(y_values, threshold=20):
    consecutives = get_consecutive_integer_series(y_values)
    if consecutives:
        distances = np.array([y2-y1  for y1,y2 in consecutives])
        if any(distance > threshold for distance in distances):
            max_index = np.argmax(distances)
            return consecutives[max_index]
        else:
            return False
    else:
        return False

def has_mid_line(image,point1,point2):
    left_point  = point1[0] - 10,round((point2[1] + point1[1])/2)
    right_point = point1[0] + 10,round((point2[1] + point1[1])/2)
    if image.getpixel(left_point) != 0 or image.getpixel(right_point) != 0:
        return True
    else:
        return False

def pixel_to_mumeter(edges_coordinates,shape,radius):
    ############# determine scale ##################
    edges_array = np.zeros(shape)
    for x,y in edges_coordinates:
        edges_array[y][x] = 255

    edges_img=Image.fromarray(edges_array)

    # Organize coordinates by x value
    coordinates_by_x = defaultdict(list)
    for x, y in edges_coordinates:
        coordinates_by_x[x].append(y)

    # Identify vertical lines
    vertical_lines = []
    for x, y_values in coordinates_by_x.items():
        consecutives = could_form_vertical_line(y_values)
        if consecutives and has_mid_line(edges_img,(x,consecutives[0]),(x,consecutives[1])):
            vertical_lines.append((x, y_values))

    # Calculate the absolute differences between x-coordinates of all pairs
    diffs = [abs(a[0] - b[0]) for a in vertical_lines for b in vertical_lines if a != b]
    # Find the maximum difference -> number of pixels corresponding to scale legend
    max_diff = max(diffs)
    scale = 1000 if max_diff > 2000 else 500
    square_radius = int((max_diff/scale) * radius)
    return square_radius

def score_function(file1,file2,folder,radius = 10):
    image1 = Image.open(folder + '/' + file1).convert('L') # convert to grayscale
    image2 = Image.open(folder + '/' + file2).convert('L')

    pixels1 = np.asarray(image1)
    pixels2 = np.asarray(image2)

    u1 = np.zeros(pixels1.shape) # initialize np.array
    u2 = np.zeros(pixels2.shape)
    edges = []
    for x in range(len(pixels1[1])):
        for y in range(len(pixels1)):
            if pixels1[y][x] != pixels2[y][x]: # filter out same pixels borders texts etc
                u1[y][x] = pixels1[y][x] # seg1 pixels go to u1
                u2[y][x] = pixels2[y][x] # seg2 pixels go to u2
            else:
                if image1.getpixel((x,y)) != 0:
                    edges.append((x,y))

    seg1_coords = np.nonzero(u1) # get coordinates where u1 is not zero
    seg2_coords = np.nonzero(u2) # get coordinates where u1 is not zero

    #labels = ['CD68'] * len(seg1_coords[0]) + ['CK'] * len(seg2_coords[0])
    seg1_coords_zipped = [[x,y] for x,y in zip(seg1_coords[0],seg1_coords[1])]
    seg2_coords_zipped = [[x,y] for x,y in zip(seg2_coords[0],seg2_coords[1])]

    square_radius = pixel_to_mumeter(edges,pixels1.shape,radius)

    positives = 0
    for x,y in seg1_coords_zipped:
        subset_seg2 = u2[x - square_radius:x + square_radius , y -square_radius:y+square_radius]
        if np.all(subset_seg2 == 0):
            continue
        else:
            positives += 1

    coords_zipped = seg1_coords_zipped + seg2_coords_zipped
    score = (positives*100)/len(seg1_coords[0])
    #score = silhouette_score(coords_zipped, labels,sample_size = 10000,metric= 'manhattan')
    #score = calinski_harabasz_score(coords_zipped, labels)
    #score = davies_bouldin_score(coords_zipped, labels)
    return score,seg1_coords,seg2_coords

def plot_roi(res):
    folder = res['folder']
    roi = res['roi']
    fig, axs = plt.subplots(1, 3, gridspec_kw={'width_ratios': [4, 1, 1]})
    axs[0].scatter(res['coords']['CD68'][1],res['coords']['CD68'][0],color = cd68_color ,s = 0.3)
    if 'CD45' in res['coords'] :
        score_text = "{}%".format(round(res['scores']['CD68']['cd45_score'],1))
        axs[0].scatter(res['coords']['CD45'][1],res['coords']['CD45'][0],color = cd45_color,s = 0.3)
        axs[2].bar(0,height = res['scores']['CD68']['cd45_score'],color = cd45_color)
        axs[2].text(x = 0,y = res['scores']['CD68']['cd45_score'] + 5 ,s = score_text, ha='center', va='top',fontsize = 18)

    if 'CK' in res['coords']:
        score_text = "{}%".format(round(res['scores']['CD68']['ck_score'],1))
        axs[0].scatter(res['coords']['CK'][1],res['coords']['CK'][0],color = ck_color,s = 0.3)
        axs[1].bar(0,height = res['scores']['CD68']['ck_score'],color = ck_color)
        axs[1].text(x = 0,y = res['scores']['CD68']['ck_score'] + 5 ,s = score_text, ha='center', va='top',fontsize = 18)

    axs[0].axis('off')

    axs[1].set_ylim(0,80)
    axs[1].axis('off')

    axs[2].set_ylim(0,80)
    axs[2].axis('off')
    # Adjust the layout and display the plot
    plt.tight_layout()
    #plt.title("roi:{}  folder:{}".format(roi,folder[:15]))
    plt.savefig(os.path.join(geomx_images,"results", "{}_{}.png".format(folder,roi)))
    plt.close()
    return None

def slide_function(folder):
    roi_pattern = r'- (0\d{2})'
    files = os.listdir(folder)
    # Dictionary to group matches
    roi_matches = {}
    # Iterate over each string in the list
    for file in files:
        # Find all matches in the current string
        matches = re.finditer(roi_pattern, file)
        # Iterate over each match
        for match in matches:
            # Extract the group (the number part)
            roi = match.group(1)
            # If the group is not in the dictionary, add it with an empty list
            if roi not in roi_matches:
                roi_matches[roi] = []
            # Append the current string to the group's list
            roi_matches[roi].append(file)

    triplettes = []
    for roi, items in roi_matches.items():
        cd68_index = [i for i in items if re.search(r'CD68',i)]
        ck_index = [i for i in items if re.search(r'CK',i)]
        cd45_index = [i for i in items if re.search(r'CD45',i)]
        triplettes.append((cd68_index,ck_index,cd45_index,roi))
    cleaned = lambda x : '' if len(x) == 0 else x[0]

    triplettes = [(cleaned(x),cleaned(y),cleaned(z),roi) for x,y,z,roi in triplettes]
    results = []
    for index, (cd68_seg,ck_seg,cd45_seg,roi) in enumerate(triplettes):
        scores = {}
        coords = {}
        if len(cd68_seg) > 0 :
            scores['CD68'] = {'ck_score':0,'cd45_score':0}
            if len(ck_seg) > 0:
                scores['CD68']['ck_score'] , coords['CD68'],coords['CK'] = score_function(cd68_seg,ck_seg,folder,radius = 10)
            if len(cd45_seg)  > 0:
                scores['CD68']['cd45_score'] , coords['CD68'],coords['CD45'] = score_function(cd68_seg,cd45_seg,folder,radius = 10)
        if len(cd45_seg) > 0 :
           scores['CD45'] = {'cd68_score':0,'ck_score':0}
           if len(cd68_seg) > 0:
               scores['CD45']['cd68_score'] , coords['CD45'],coords['CD68'] = score_function(cd45_seg,cd68_seg,folder,radius = 10)
           if len(ck_seg)  > 0:
               scores['CD45']['ck_score'] , coords['CD45'],coords['CK'] = score_function(cd45_seg,ck_seg,folder,radius = 10)

        if len(ck_seg) > 0 :
           scores['CK'] = {'cd68_score':0,'cd45_score':0}
           if len(cd68_seg) > 0:
               scores['CK']['cd68_score'] , coords['CK'],coords['CD68'] = score_function(ck_seg,cd68_seg,folder,radius = 10)
           if len(cd45_seg)  > 0:
               scores['CK']['cd45_score'] , coords['CK'],coords['CD45'] = score_function(ck_seg,cd45_seg,folder,radius = 10)

        res = {
               'scores':scores,
               'coords':coords,
               'roi'   :roi,
               'folder':os.path.basename(folder)
               }
        results.append(res)
    return results

annotations = pd.read_csv(annotations_file)

folder_slide_dict = {
        "18H01000.12_20H10335.6_SEIN_L02_S01":"18H01000.12_20H10335.6_SEIN_L02",
        "18H03710-3 primitive colon _3 11-10-22":"18H03710-3 primitive colon #3 11-10-22",
        "19H07880.10_21H09355.7_SEIN_L02_S01":"19H07880.10_21H09355.7_SEIN_L02",
        "20H08710-11 primitive colon _1 13-10-22":"20H08710-11 primitive colon #1 13-10-22",
        "20H10907.05 PRIMITIVE 07-07-2022":"20H10907.05 PRIMITIVE 07-07-2022",
        "21H11316.08_21H00501.09_COLON_L02_S01":"21H11316.08_21H00501.09_COLON_L02",
        "22H01374.14_21H05510.11_COLON_L02_S01":"22H01374.14_21H05510.11_COLON_L02",
        "22H01892-17 primitive colon_4 06-10-22":"22H01892-17 primitive colon#4 06-10-22"
        }

print(1)
with Pool(processes = 8) as pool:
    results = pool.map(
                slide_function,
                [geomx_images + '/' + x for x in list(folder_slide_dict.keys())])
print(2)
results = [item for sublist in results for item in sublist]
print(3)
results = sorted(results,key = sort_key)
with open('results.pickle', 'wb') as file:
        pickle.dump(results, file)
print(results[0])

print(4)
for res in [d for d in results if 'CD68' in d.get('coords', {})]:
    plot_roi(res)

print(5)
images = [x for x in glob.glob(geomx_images+'/results/*') if x.endswith(".png")]
combine_images(columns=8, space=0, images=images)
print(6)

############ update annotations ##############
annotations['ck_score'] = None
annotations['cd45_score'] = None
annotations['cd68_score'] = None
for i in results:
    slide = folder_slide_dict[i['folder']]
    roi = i['roi']
    segments = list(i['scores'].keys())
    for segment in segments:
        mask = (annotations['slide name'] == slide) & (annotations['roi'] == float(roi)) & (annotations['segment'] == segment)
        if segment == 'CD68':
            annotations.loc[mask, 'cd45_score'] = i['scores'][segment]['cd45_score']
            annotations.loc[mask, 'ck_score'] = i['scores'][segment]['ck_score']
        if segment == 'CK':
            annotations.loc[mask, 'cd45_score'] = i['scores'][segment]['cd45_score']
            annotations.loc[mask, 'cd68_score'] = i['scores'][segment]['cd68_score']
        if segment == 'CD45':
            annotations.loc[mask, 'ck_score'] = i['scores'][segment]['ck_score']
            annotations.loc[mask, 'cd68_score'] = i['scores'][segment]['cd68_score']

print(annotations.to_string())
annotations.to_csv("/home/dimitris/Workspace/crc/data/geomx_annotations_with_score.csv")
