# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 14:03:56 2018

@author: vlieg
"""


########################

## This script runs all the neccessary functions defined in 
## "functions_2Pspectra.py" to fit the excitation spectra of emitters.
##
## Step 1: Find emitters within the 2D image and save coordinates (global_coordinates)
## Step 2: Get sub-pixel location accuracy by fitting a Gaussian to the 
##         found global coordinates (local_coordinates) and a Squared Lorentzian
##         to acquire the spectral parameters 
## Step 4: Concatonate all global coordinates and fit parameters in a DataFrame 
##         and write to a .xlsx file. 

########################

#%% Load Necessary Libraries
import numpy as np
import file_io as file_io
import functions_2Pspectra as func
import pandas as pd
import matplotlib.pyplot as plt

#%% Get filename & File parameters 
file_path = file_io.get_path()
file_log  = file_io.read_logfile(file_path)

ypix   = np.uint(str.split(file_log[8],'=')[1])
xpix   = np.uint(str.split(file_log[9],'=')[1])
zsteps = file_log[18]
zsteps = str.split(zsteps,"=")[1]
zsteps = np.uint(str.split(zsteps,",")[0])

if zsteps == 1:
    nstacks=1
else:
    nframes = np.uint(str.split(file_log[7],'=')[1])
    nstacks = int(nframes/zsteps)


#%% Get Global coordinates from stack  
threshold_factor = 1.3
mask_size = [11,11,zsteps]
max_num_peaks = 0

global_coords,num_traces_pstack = func.get_global_coords(file_path,threshold_factor,mask_size,max_num_peaks)


#%% 3D Fit (x,y,lambda) to found Volume-Of-Interest by Global Coordinates   
for stack_nr in range(0,nstacks):
    pd_fitparams=func.get_local_and_spectrum(file_path,global_coords,mask_size,stack_nr)


#%% Write Fit Parameters to .xlsx file 
file_io.write_pd2xlsx(file_path,pd_fitparams)
