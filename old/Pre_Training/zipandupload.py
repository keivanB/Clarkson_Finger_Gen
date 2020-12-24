#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 16:03:54 2019

@author: This code zips  Precise dataset for chamelon
"""
import os
import zipfile
import numpy as np
import glob

inputloc = '/media/citer/My Passport/Precise Dataset/Crossmatch Guardian'
folder_list=sorted(os.listdir(inputloc))
folder_list.remove('Dataset Description.pdf')
outputloc = '/media/citer/My Passport/Precise Dataset/zips'

for folder in folder_list:
    inputf = os.path.join(inputloc, folder, 'raw', 'live')
    file_list = glob.glob(inputf+'/*/*')
    fzip = zipfile.ZipFile(os.path.join(outputloc, folder+'.zip'), 'w')
    for file in file_list:
        fzip.write(file, file, compress_type = zipfile.ZIP_DEFLATED)
    fzip.close()
