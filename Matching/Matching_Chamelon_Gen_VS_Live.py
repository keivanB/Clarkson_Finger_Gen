import os
import subprocess
import multiprocessing
import concurrent.futures
import numpy as np
      
def Matching(fn):
    p = subprocess.check_call(["/NBIS/bozorth3/bin/bozorth3", "-o", "/Synt/Matching_Real/"+str(fn)+".scr", "-A", "outfmt=sg", "-A", "maxfiles=73000", "-p", "/Synt/Gen_Minuteas/"+str(fn)+".xyt", "-G", "/Synt/Real_Live.txt"]) 
    
file_dic = {}
nprocess = 64
Nfold =  int(50048/nprocess)

Nnodes = 4
C = int(Nfold/Nnodes)

node = 4

if node !=4:
    print('calculaiting folds from {} to {} out of {}'.format((node-1)*C, (node)*C, Nfold))
    for fold in range((node-1)*C, (node)*C):
         file_dic[fold]=range((fold)*nprocess,(fold+1)*nprocess)

    for key in file_dic.keys():
        with concurrent.futures.ProcessPoolExecutor(nprocess) as executor:
            results = [executor.submit(Matching, file) for file in file_dic[key]]

if node ==4:
    print('calculaiting folds from {} to {} out of {}'.format((node-1)*C, Nfold, Nfold))
    for fold in range((node-1)*C, Nfold):
         file_dic[fold]=range((fold)*nprocess,(fold+1)*nprocess)

    for key in file_dic.keys():
        with concurrent.futures.ProcessPoolExecutor(nprocess) as executor:
            results = [executor.submit(Matching, file) for file in file_dic[key]]