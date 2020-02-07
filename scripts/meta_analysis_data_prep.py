import pandas as pd
import os
import re
from functools import reduce

file_list = os.listdir('/data/vitaled2/summary_stats')
suffixes = [re.search(r'\.(.*?)\.', cohort).group(1) for cohort in file_list]

summary_data_list = [pd.read_csv('/data/vitaled2/summary_stats/' + file, sep='\t') for file in file_list]

for i in range(len(summary_data_list)):
    summary_data_list[i]['cohort'] = suffixes[i]

summary_df = pd.concat(summary_data_list)

# drop duplicate markers within cohorts, keeping maker with lowest pval
summary_df_cleaned = summary_df.sort_values('P').drop_duplicates(['markerID','cohort'], keep='first')



# markers_list = []
# for i in range(len(file_list)):
#     markers = pd.read_csv('summary_stats/' + file_list[i], sep='\t', usecols=['markerID'])
#     markers_list.append(summary)

# select markers that exist across all datasets
merged = reduce(lambda x, y: pd.merge(x, y, on=['markerID'], how='inner'), summary_data_list) 
marker_list = list(merged.drop(merged.columns.difference(['markerID']), axis=1).drop_duplicates('markerID').markerID)

final_df = summary_df_cleaned[summary_df_cleaned.markerID.isin(marker_list)]
final_sorted = final_df.sort_values(['markerID','cohort'])
final = final_sorted[['markerID', 'cohort', 'beta', 'se']]

final.to_csv("random_effect_test_5_cohorts.csv", index=False, sep='\t')
