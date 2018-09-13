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
##         found global coordinates (local_coordinates)
## Step 3: Get spectra parameters by fitting a squared Lorentzian to the found 
##         global coordinates.
## Step 4: Concatonate all global coordinates and fit parameters in a DataFrame 
##         and write to a .xlsx file. 

########################

#%% Load Necessary Libraries
import numpy as np
import file_io as file_io
import functions_2Pspectra as func 


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

#Remove traces with wrong coordinates
if len(global_coords) == 1:
    global_coords=global_coords[0]
    
global_coords_corrected = np.empty([0,3])
trace_nr=0
for trace_coords in global_coords:
    if trace_coords[0]<5 or trace_coords[1]<5:
        trace_nr=trace_nr+1
    else:
        trace_coords = np.expand_dims(trace_coords,1)
        trace_coords = np.swapaxes(trace_coords,0,1)
        global_coords_corrected=np.append(global_coords_corrected,trace_coords,0)
        trace_nr=trace_nr+1
        
for stack_nr in range(0,nstacks):
    local_coordinates,fit_parameters_coordinates,fit_errors_coordinates = func.get_local_coordinates(file_path,global_coords_corrected,stack_nr,mask_size)



#%% Fit Lorentzian to Excitation Spectra Data 
fit_parameters_spectra, fit_errors_spectra = func.fit_spectra(file_path,global_coords_corrected,mask_size)


#%% Combine Local Coordinates and Spectra fit Parameters and write to .xlsx file
#Insert the local coordinates in the array of the 3D Gaussian Fit Parameters  
fit_parameters_coordinates_wlocal=fit_parameters_coordinates
fit_parameters_coordinates_wlocal[:,2:5]=local_coordinates

traces_parameters = np.array(np.concatenate((global_coords_corrected,fit_parameters_coordinates,fit_parameters_spectra,fit_errors_coordinates,fit_errors_spectra),1))
column_headers= ['x_global (pix)',
                 'y_global (pix)',
                 'z_global (pix)',
                 'amplitude_2D (a.u.)',
                 'offset_2D (a.u.)',
                 'x_local (pix)',
                 'y_local (pix)',
                 'z_local (pix)',
                 'width_x (pix)',
                 'width_y (pix)', 
                 'width_z (pix)', 
                 'spectrum peak (nm)',
                 'spectrum width (nm)',
                 'spectrum height (a.u.)',
                 'offset_spectrum (a.u.)',
                 'err_amplitude_2D (a.u.)',
                 'err_offset_2D (a.u.)',
                 'err x_local (pix)',
                 'err_y_local (pix)',
                 'err_z_local (pix)',
                 'err_width_x (pix)', 
                 'err_width_y (pix)', 
                 'err_width_z (pix)',
                 'err_spectrum peak (nm)',
                 'err_peak width (nm)',
                 'err_peak height (a.u.)',
                 'err_offset (a.u.)']

column_indices = range(0,len(column_headers))

## Write to Excel file and convert to Pandas Dataframe 
fit_parameters = file_io.write_xlsx_array(file_path,traces_parameters,column_headers,column_indices)


#%% 
