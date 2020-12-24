import os
import pandas as pd
from glob import glob
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


train_df = pd.read_pickle('/Synt/train_df_min.p')
train_df.FilePath = train_df.FilePath.apply(lambda x: x[:-3]+'scr')
input_matches_folder = '/Synt/Matching_Real_VS_Real/'
files = os.listdir(input_matches_folder)

train_samples = train_df.FilePath.tolist()
error = []
indexes = []
for i, file in enumerate(files):
    if file not in train_samples:
        error.append(file)
        indexes.append(i)

renamed = []
for j, file in enumerate(error):
    # correct files with bad naming
    if file[-3:] != 'scr':
        os.rename(os.path.join(input_matches_folder, files[indexes[j]]), os.path.join(input_matches_folder, files[indexes[j]]+'.scr'))
        renamed.append([file,j])
