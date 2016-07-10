# -*- coding: iso-8859-1 -*-
# Necessary imports 
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 


#
data_source = 'IN2P3'
dir_info_files = '/astro/store/scratch/tmp/suberlak/S13Agg/repo_fls/'
name_detection = 'DeepSource'+data_source+'_i_lt235.csv.gz'
detection_data = pd.read_csv(dir_info_files+name_detection)
narrow_detection_data  = detection_data[['deepSourceId', 'ra', 'decl']]
narrow_detection_data.to_csv(dir_info_files+ name_detection[:-7]+'_very_narrow.csv')