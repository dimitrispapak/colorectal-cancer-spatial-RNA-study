import pandas as pd
import numpy as np
excel_labworksheet = pd.read_excel('/mnt/beegfs/userdata/d_papakonstantinou/crc/geomx_crc/ROI_GeoMx_RNA_WTA_NGS_RT16.22.xlsx',sheet_name='LabWorsheet RT016.22',skiprows = 9, dtype={'Roi': str})
excel_labworksheet = excel_labworksheet.rename(columns = {'Slide Name':'slide name','Scan Name':'scan name','Panel':'panel','Roi':'roi','Segment':'segment','Aoi':'aoi','Area':'area'})
excel_labworksheet['mucinous'] = 'true'
excel_labworksheet_subset = excel_labworksheet[['Sample_ID', 'slide name', 'panel', 'roi', 'segment','aoi', 'area','mucinous']]

old_annotations = pd.read_excel('/mnt/beegfs/userdata/d_papakonstantinou/geomx/output.xlsx', dtype={'roi': str})

old_annotations = old_annotations[old_annotations.pancreas_colon != 'pancreas']
old_annotations['tissue'] = np.where(old_annotations['primitive_metastasis']== 'primitive','colon',old_annotations['metastasis_location'])
old_annotations['mucinous'] = 'false'
old_annotations_subset = old_annotations[['Sample_ID', 'slide name', 'panel', 'roi', 'segment','aoi', 'area', 'tissue', 'patient_id', 'location','mucinous']]

excel_metadata = pd.read_excel('/mnt/beegfs/userdata/d_papakonstantinou/crc/geomx_crc/ROI_GeoMx_RNA_WTA_NGS_RT16.22.xlsx',sheet_name='ROI RT016.22',skiprows = 37)
excel_metadata = excel_metadata.rename(columns = {
	"Echantillon":"patient_id",
	"Nom lame \nsur GeoMx":"slide name",
	"Groupe":"tissue",
	"N° ROI":"roi",
	"Nom du ROI":"location",
	'AOI':"aoi",
	'Nom AOI':"aoi type",
	'Marqueurs \nsegmentés':"segment",
	'Aire AOI\nréelle \npost-run (µm2)':"area",
	'Nombre de \nNoyaux AOI':"nuclei"
	})
excel_metadata = excel_metadata[["slide name","patient_id","tissue","roi","location","aoi","aoi type","segment","area","nuclei"]]
excel_metadata = excel_metadata[excel_metadata['slide name'].isin(excel_labworksheet_subset['slide name'].to_list())]
fix_roi_name = lambda x: "0"+ x.split('_')[1]
excel_metadata['roi'] = excel_metadata['roi'].apply(fix_roi_name)

def rename_location(value):
	if value == 'Tumor':
		return 'intra-tumor'
	elif value == 'Normal':
		return 'extra-tumor'
	else:
		return 'front'

excel_metadata['location'] = excel_metadata['location'].apply(lambda x: rename_location(x))
excel_metadata['tissue'] = excel_metadata['tissue'].apply(lambda x: x.lower())

result = pd.concat([excel_labworksheet_subset, old_annotations_subset], ignore_index=True, sort=False)
for index , row in result.iterrows():
	slide = row['slide name']
	roi = row['roi']
	segment = row['segment']
	subset = excel_metadata[(excel_metadata['slide name'] == slide) & \
				(excel_metadata['roi'] == roi) & \
				(excel_metadata['segment'] == segment)]

	if len(subset) == 1 :
		dict_ = subset.iloc[0].to_dict()
		result.loc[index,'patient_id'] = dict_['patient_id']
		result.loc[index,'tissue'] = dict_['tissue']
		result.loc[index,'location'] = dict_['location']


names = {'colon_patient_1':'patient 1',
'colon_patient_2':'patient 2',
'colon_patient_3':'patient 3',
'colon_patient_4':'patient 4',
'21H11316.08':'patient 5',
'21H00501.09':'patient 6',
'22H01374.14':'patient 7',
'21H05510.11':'patient 8',
'18H01000.12':'patient 9',
'21H09335.07':'patient 10',
'19H07880.10':'patient 11',
'20H10335.6':'patient 12'}

result['patient_id'] = result['patient_id'].map(names)
# only uncomment if you know what you are doing!
result.to_excel('/mnt/beegfs/userdata/d_papakonstantinou/crc/geomx_crc/annotations.xlsx',index = False)
