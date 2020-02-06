import pandas as pd
import os
import re
from functools import reduce

file_list = os.listdir('summary_stats')
suffixes = [re.search(r'\.(.*?)\.', cohort).group(1) for cohort in os.listdir('summary_stats')]

# summary_data_list = [pd.read_csv(file, sep="\t") for file in file_list]

summary_list = []
for i in range(len(file_list)):
    summary = pd.read_csv('summary_stats/' + file_list[i], sep='\t', usecols=['markerID', 'minorAllele', 'majorAllele', 'beta', 'se'])
    summary.columns = ['markerID', 'minorAllele', 'majorAllele', 'beta_' + suffixes[i], 'se_' + suffixes[i]]
    summary_list.append(summary)
    
merged = reduce(lambda x, y: pd.merge(x, y, on=['markerID', 'minorAllele', 'majorAllele'], how='inner'), summary_list)
  