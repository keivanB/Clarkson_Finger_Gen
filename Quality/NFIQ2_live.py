#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 18:53:14 2019

@author: Keivan Bahmani

        This Code Extract the NFIQ@ Quality scores using NFIQ2K1 docker image
        
        Keivan Bahmani
"""
# =============================================================================
# Libraries
# =============================================================================
import os
#from PIL import Image
#import cv2
#import shutil
# =============================================================================
Docker = True
if Docker ==True:
    input_folder = '/Synt/Live_512_BMP/gray'
    output_folder = '/Synt/NFIQ2_Live'

wrongfile=[]
if not os.path.exists(os.path.join(output_folder)):
    os.makedirs(os.path.join(output_folder))
List_files = os.listdir(os.path.join(input_folder))
os.chdir(os.path.join(output_folder))

for file in List_files:
    command ='/NFIQ2/NFIQ2/bin/NFIQ2 SINGLE '+os.path.join(input_folder, file)+' BMP true true >'+os.path.join(output_folder, file[:-3]+'txt')
    w = os.system(command)
    if w!=0:
        wrongfile.append(os.path.join(file))
    print(file)
file=open(os.path.join(output_folder,'wrongfile.txt'), 'w')
for line in wrongfile:
    file.write("%s\n" %line)