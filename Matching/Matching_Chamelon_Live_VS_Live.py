import os
import subprocess
import multiprocessing
import concurrent.futures
import numpy as np
import pandas as pd
      
def Matching(fn):
    p = subprocess.check_call(["/NBIS/bozorth3/bin/bozorth3", "-o", "/Synt/Matching_Real_VS_Real/"+str(fn)+".scr", "-A", "outfmt=sg", "-A", "maxfiles=73000", "-p", "/Synt/Live_Minuteas/"+str(fn)+".xyt", "-G", "/Synt/Real_Live.txt"]) 
    
file_dic = {}
nprocess = 64


Train_df = pd.read_pickle('/Synt/train_df_min.p')
Train_df.FilePath = Train_df.FilePath.apply(lambda x: x[0:-4])
files = Train_df.FilePath.to_list()

Nfold =  int(len(files)/nprocess)+1

Nnodes = 4
             
C = int(Nfold/Nnodes)

node = 4

if node !=4:
    print('calculaiting folds from {} to {} out of {}'.format((node-1)*C, (node)*C, Nfold))
    for fold in range((node-1)*C, (node)*C):
         file_dic[fold]=files[(fold)*nprocess:(fold+1)*nprocess]

    for key in file_dic.keys():
        with concurrent.futures.ProcessPoolExecutor(nprocess) as executor:
            results = [executor.submit(Matching, file) for file in file_dic[key]]

if node == 4:
    print('calculaiting folds from {} to {} out of {}'.format((node-1)*C, Nfold, Nfold))
    for fold in range((node-1)*C, Nfold-1):
         file_dic[fold]= files[(fold)*nprocess:(fold+1)*nprocess]        
    file_dic[Nfold-1]= files[(Nfold-1)*nprocess:len(files)]     

    for key in file_dic.keys():
        with concurrent.futures.ProcessPoolExecutor(nprocess) as executor:
            results = [executor.submit(Matching, file) for file in file_dic[key]]