import pandas as pd
import os
import shutil
import os.path

OLD='/mnt/beegfs/userdata/d_papakonstantinou/geomx/output.xlsx'
OLD_DCCS='/mnt/beegfs/userdata/d_papakonstantinou/geomx/dccs'
NEW='/mnt/beegfs/userdata/d_papakonstantinou/crc/geomx_crc/RT016.22_RNA_WTA_20240409T1217/RT016.22_RNA_WTA_20240409T1217_LabWorksheet.txt'
new = pd.read_csv(NEW,sep = '\t',skiprows =16)
old = pd.read_excel(OLD)
sample_ids = old.Sample_ID[old['pancreas_colon'] != 'pancreas'].to_list()

for id  in  sample_ids :
	src = '/mnt/beegfs/userdata/d_papakonstantinou/geomx/dccs/' + id + '.dcc'
	dst = '/mnt/beegfs/userdata/d_papakonstantinou/crc/geomx_crc/dccs_2/'
	shutil.copy(src, dst)
