#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 00:09:25 2019

@author: Keivan Bahmani

        Crop the ITR Crossmatch section to 480by 480 for the GAN!!!
"""
import os 
import cv2 as cv
import matplotlib.pyplot as plt
inputfolder = '/media/citer/Data/Crossmatch/Live/CrossMatch'
outputfolder = '/media/citer/Data/Chrosmatch/Resized'

if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)

folders = os.listdir(inputfolder)
folders = ['9']
for f in folders:
    inp = os.path.join(inputfolder,f)
    File_list = os.listdir(inp)
    if 'Thumbs.db' in File_list:
        File_list.remove('Thumbs.db')
    for file in File_list:
        img  = cv.imread(os.path.join(inp,file))
        cv.imwrite(os.path.join(outputfolder, file),img[:,80:560])
