import pickle
import pandas as pd
import numpy as np

#pd.options.mode.chained_assignment = None

annotations_file = '/mnt/backup/Workspace/crc/data/geomx_annotations.csv'
with (open("results.pickle", "rb")) as openfile:
	results = pickle.load(openfile)

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

annotations = pd.read_csv(annotations_file)
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



none_rows = annotations[annotations[['ck_score', 'cd45_score', 'cd68_score']].isna().all(axis=1)]
full_rows = annotations[annotations[['ck_score', 'cd45_score', 'cd68_score']].notna().any(axis=1)]
print(len(annotations))
print(len(none_rows))
print(len(full_rows))
print(*sorted(set(none_rows['slide name'])),sep= "\n")
print(*sorted(set(full_rows['slide name'])),sep= "\n")
print(set(none_rows['slide name']).intersection(set(full_rows['slide name'])))
