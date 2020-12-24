import os
import subprocess
import multiprocessing
import concurrent.futures
import numpy as np
      
def Matching(fn):
    p = subprocess.check_call(["/NBIS/bozorth3/bin/bozorth3", "-o", "/Synt/Matching_Gen/"+str(fn)+".scr", "-A", "outfmt=sg", "-A", "maxfiles=51000", "-p", "/Synt/Gen_Minuteas/"+str(fn)+".xyt", "-G", "/Synt/matching_lists/All_Gen.txt"]) 
    
file_dic = {}
nprocess = 8
Nfold =  int(50048/nprocess)

for fold in range(Nfold):
     file_dic[fold]=range((fold)*nprocess,(fold+1)*nprocess)

for key in file_dic.keys():
    with concurrent.futures.ProcessPoolExecutor(nprocess) as executor:
        results = [executor.submit(Matching, file) for file in file_dic[key]]
