import os
import subprocess
import multiprocessing

      
def Start_Sampling(fn):
    p = subprocess.Popen(["/NBIS/bozorth3/bin/bozorth3", "-o", "/Synt/matching/"+str(fn)+".scr", "-A", "outfmt=sg", "-p", "/Synt/Gen_Minuteas/1.xyt", "-G", "/Synt/matching_lists/0.txt"])        

   
Process_List = []
Inputs = []
samples = 10
for i in range(samples):
    Inputs.append([int(i)])

for i in range(n_process):
    P = multiprocessing.Process(target=Start_Sampling, args=[Inputs[i][0]])
    P.start()
    Process_List.append(P)

for P in Process_List:
    P.join()
