#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 23:18:30 2019

@author: Keivan Bahmani

        PCA of the Fingerprint !!!
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import cv2 as cv
import os 

size = 480*480
inputfolder = '/media/citer/Data/Crossmatch/Images/Resized'
File_List = os.listdir(inputfolder)
file_list = File_List[0:1000]
Images =[]
for file in file_list :    
    img = cv.imread(os.path.join(inputfolder, file),0)
    Images.append(np.reshape(img, size))
Images = np.array(Images)
Images_sc = StandardScaler().fit_transform(Images)

pca = PCA(512)
Images_PCA = pca.fit_transform(Images_sc)
Images_final = pca.inverse_transform(Images_PCA)

np.sum(pca.explained_variance_ratio_)  # going from 512x512 images to Latent space of 512 we can explain 92% of the varience

c1 = np.reshape(Images_final[1,:], (800,600))
c2 = np.reshape(Images_final[2,:], (800,600))
c3 = np.reshape(Images_final[3,:], (800,600))
plt.figure()
plt.imshow(c1, cmap='gray')
plt.figure()
plt.imshow(np.reshape(Images[1,:], (800,600)), cmap='gray')

plt.imshow(c2, cmap='gray')
plt.figure()
plt.imshow(c3, cmap='gray')

