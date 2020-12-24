#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 18:53:14 2019

@author: Keivan Bahmani

        This Code extract the minutea points from the fingerprints using the
        NBISK2 docker image
        
        Keivan Bahmani
"""
# =============================================================================
# Libraries
# =============================================================================
import os
#from PIL import Image
#import cv2
#import shutil

Docker = True
if Docker ==True:
    input_folder = '/Synt/gray'
    output_folder = '/Synt/Minuteas'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
List_files = os.listdir(input_folder)
os.chdir(output_folder)

for file in List_files:
    if file[-4:]=='JPEG':
        command='/nbis/bin/mindtct'+' '+os.path.join(input_folder,file)+' '+' '+file[:-5]
        os.system(command)
        print(file)
