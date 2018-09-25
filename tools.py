# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:25:04 2018

@author: vlieg
"""



# =============================================================================
# Tools & Legacy VIs
# 
# 
# =============================================================================

import numpy as np
import tkinter as tk
import file_io as file_io
from scipy.optimize import curve_fit
from tkinter import filedialog
import pandas as pd 
import functions_2Pspectra as func

#%% Get max value in array with corresponding indices
def get_max(array):
    indices = np.unravel_index(np.argmax(array), np.shape(array))
    max_intensity = array[indices]
    indices = np.int16(indices)    
    return max_intensity, indices
    
    
#%% Get min value in array with corresponding indices
def get_min(array):
    indices = np.unravel_index(np.argmin(array), np.shape(array))
    max_intensity = array[indices]
    indices = np.int16(indices)
    return max_intensity, indices


#%% Get Peaks from Stack 
def create_3D_mask(size_x,size_y,size_z,amp,sigma):
    # Create 3D Gaussian mask 
    
    A = amp   # amplitude
    S = sigma # sigma
    
    x = np.arange(0,size_x,1)
    y = np.arange(0,size_y,1)
    z = np.arange(0,size_z,1)
    mid_indx = np.uint16(np.floor(len(x)/2))
    mid_indy = np.uint16(np.floor(len(y)/2))
    mid_indz = np.uint16(np.floor(len(z)/2))
    
    # Create xyz-coordinates 
    xx, yy, zz = np.meshgrid(x, y, z, sparse=False)
    
    # Calculate 3D Gaussian
    gauss_x = np.square(xx-x[mid_indx])/(2*np.square(S))
    gauss_y = np.square(yy-y[mid_indy])/(2*np.square(S))
    gauss_z = np.square(zz-z[mid_indz])/(2*np.square(S))
    gauss_xyz = A*np.exp(-1*(gauss_x+gauss_y+gauss_z))
    
    # Get Mask in boolean
    mask_xyz = gauss_xyz<np.median(gauss_xyz)
    
    return mask_xyz

#%% Create a 2D mask based on a 2D Gaussian
def create_2D_mask(size_x,size_y,amp,sigma):
    # Create 3D Gaussian mask 
    
    A = amp   # amplitude
    S = sigma # sigma
    
    x = np.arange(0,size_x,1)
    y = np.arange(0,size_y,1)
    mid_indx = np.uint16(np.floor(len(x)/2))
    mid_indy = np.uint16(np.floor(len(y)/2))
    
    # Create xyz-coordinates 
    xx, yy = np.meshgrid(x, y, sparse=False)
    
    # Calculate 3D Gaussian
    gauss_x = np.square(xx-x[mid_indx])/(2*np.square(S))
    gauss_y = np.square(yy-y[mid_indy])/(2*np.square(S))
    gauss_xy = A*np.exp(-1*(gauss_x+gauss_y))
    
    # Get Mask in boolean
    mask_xy = gauss_xy<np.median(gauss_xy)
    
    return mask_xy

#%% Model function to be used to fit to the data:
    
def gauss_3D(xyz, p):
    xx,yy,zz = xyz[0],xyz[1],xyz[2]
    A,C,x0,y0,z0,wx,wy,wz = p
    gauss_x = np.square(xx-x0)/(2*np.square(wx))
    gauss_y = np.square(yy-y0)/(2*np.square(wy))
    gauss_z = np.square(zz-z0)/(2*np.square(wz))
    
    return A*np.exp(-1*(gauss_x+gauss_y+gauss_z))+C


#%% Model function to be used to fit to the data:
    
def gauss_2D(xy, *p):
    xx,yy = xy[0],xy[1]
    A,C,x0,y0,wx,wy = p
    gauss_x = np.square(xx-x0)/(2*np.square(wx))
    gauss_y = np.square(yy-y0)/(2*np.square(wy))
    
    return A*np.exp(-1*(gauss_x+gauss_y))+C

#%% Model function to be used to fit to the data:

def sq_lorentzian(l, *p):
    l0,gamma,C,A = p
    lorentz = 1/(1+4*(np.sqrt(2)-1)*np.square((l-l0)/gamma))
    
    return A*lorentz+C


#%% Get z-profile from 2D ROI 
def get_zprofile(file_path,trace_coordinate,stack,mask_size):
    
    #Create mask
    size_x,size_y=mask_size[0:2]
    mask = create_2D_mask(size_x,size_y,10,10)

    # Get log-file parameters     
    logfile = file_io.read_logfile(file_path)
    
    ypix = np.uint(str.split(logfile[8],'=')[1])
    xpix = np.uint(str.split(logfile[9],'=')[1])
    size_x,size_y=mask.shape
    
        
    ### Get z-profile for fitting Lorentzian squared     
    x_range = np.int16([trace_coordinate[0]-np.floor(size_x/2),trace_coordinate[0]+np.floor(size_x/2)])
    y_range = np.int16([trace_coordinate[1]-np.floor(size_y/2),trace_coordinate[1]+np.floor(size_y/2)])
    
    mask_bounds = np.int16([np.floor(size_x/2),np.floor(size_y/2)])    
    mask_temp = mask        
    
    # X-coordinates too small:
    if  trace_coordinate[0] < mask_bounds[0]:
        x_range = np.int16([0,trace_coordinate[0]+np.floor(size_x/2)])
        mask_temp = mask_temp[mask_bounds[0]-trace_coordinate[0]::,:,:]
     
    # X-coordinates too large:
    if trace_coordinate[0] + mask_bounds[0] >= xpix:
        x_range = np.int16([trace_coordinate[0]-np.floor(size_x/2),xpix])
        ind_cut = ((x_range[0])+size_x)-xpix
        mask_temp = mask_temp[0:(size_x-ind_cut),:,:]
    
    # Y-coordinates too small:
    if trace_coordinate[1] < mask_bounds[1]:
        y_range = np.int16([0,trace_coordinate[1]+np.floor(size_y/2)])
        mask_temp = mask_temp[:,mask_bounds[1]-trace_coordinate[1]::,:]
        
    # Y-coordinates too large:             
    if trace_coordinate[1] + mask_bounds[1] >= ypix:
        y_range = np.int16([trace_coordinate[1]-np.floor(size_y/2),ypix])
        ind_cut = ((y_range[0])+size_y)-ypix
        mask_temp = mask_temp[:,0:(size_y-ind_cut),:]
        
    ROI_stack = stack[x_range[0]:x_range[1]+1,y_range[0]:y_range[1]+1,:]
    
    # Get z-profile   
    data_z = np.mean(ROI_stack,(0,1))-np.mean(stack,(0,1))

    return data_z

#%% Fit Squared Lorentzian to 1D data  
def fit_lorentzian_sq(x,data_z):
    
    #Initial Fit Parameters 
    in_l0    = np.squeeze(x[get_max(data_z)[1]])
    in_gamma = np.float64(35)
    in_C     = min(data_z)
    in_A     = max(data_z)
    
    # p0 is the initial guess for the fit coefficients
    p0 = [in_l0,in_gamma,in_C,in_A]
    
    try:
        coeff, var_matrix = curve_fit(sq_lorentzian, x, data_z, p0=p0)
        perr = np.sqrt(np.diag(var_matrix))

    except RuntimeError:
        print('Error - curve_fit failed')
        coeff = np.empty(len(p0))                
        perr  = np.empty(len(p0))
        
    fit_parameters = coeff
    fit_errors    = perr
    
    return fit_parameters, fit_errors

#%% Fit Squared Lorentzian to all found global traces in a stack
    
def fit_spectra(file_path,global_coords,mask_size):
    
    #Read in all data 
    stack=file_io.read_bin_all(file_path)  
    file_path_dat = file_io.change_extension(file_path,'.dat')
    data_dat = pd.read_csv(file_path_dat,delimiter='\t',decimal =',')
 
    it=0
    
    # Allocate memory for fit parameters and errors 
    fit_parameters = np.empty([len(global_coords),4])
    fit_errors     = np.empty([len(global_coords),4])
        
    for trace_coords in global_coords:
        print('Fit Trace #',it)
        
        # Get the raw spectrum data per trace 
        data_z = get_zprofile(file_path,trace_coords,stack,mask_size)
                
        ######  Fit data_z to a squared lorentzian #######                
        x = np.array(data_dat['Spectrum peak (nm)'])
        fit_parameters_trace, fit_errors_trace = fit_lorentzian_sq(x,data_z)
        fit_parameters[it,:]=fit_parameters_trace
        fit_errors[it,:]=fit_errors_trace
        
        it=it+1
        
    return fit_parameters, fit_errors


#%% Get sub-pixel resolution of the coordinates (local_coordiantes) from global coordinates 

def get_local_coordinates(file_path,global_coords,stack_nr,mask_size):
    # Get Measurement Parameters 
    file_log = file_io.read_logfile(file_path)
    zsteps = file_log[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])
    
    ### Get ROI from the stack using Global Coordinates and Fit 3D Gaussian

    # Make 3D Mask
    size_x,size_y,size_z=mask_size
    mask = create_3D_mask(size_x,size_y,size_z,10,10)
    
    # Define Bounds
    max_bounds=(65536,65536,size_x,size_y,size_z,size_x-1,size_y-1,size_z-1)
    #min_bounds=(0,0,0,0,0,0,0,0)
    #bounds = (min_bounds,max_bounds)
    
    ### Fit ROI to a 3D Gauss 
      
    # Read one stack from .bin file 
    if zsteps == 1:
        nstacks=1
        stack = file_io.read_bin_all(file_path)
    else:
        nframes = np.uint(str.split(file_log[7],'=')[1])
        nstacks = int(nframes/zsteps)
        stack = file_io.read_stack(file_path,stack_nr)
    
    # Allocate memory for fit parameters & errors for 1 stack 
    fit_params = np.empty([len(global_coords[:,0]),len(max_bounds)])
    fit_errors = np.empty([len(global_coords[:,0]),len(max_bounds)])
    mask_size_traces = np.empty([len(global_coords[:,0]),len(mask_size)])
    
    for trace_nr in range(0,len(global_coords)):
        
        # Get ROI
        peak_coords = global_coords[trace_nr,:]
        ROI_stack,mask_size_trace = func.get_ROI_from_stack(file_path,stack,peak_coords,mask)
        
        ##### Fit Data #####
        # Do not fit Z-coordinates when measurement is in 2D
        if zsteps ==1:
            coeff,perr = fit_2D_gauss(ROI_stack,mask_size_trace)
        else:
            coeff,perr = fit_3D_gauss(ROI_stack,mask_size_trace)
            
        fit_params[trace_nr,:]= coeff
        fit_errors[trace_nr,:]= perr
        mask_size_traces[trace_nr,:] = mask_size_trace
        
        
    ### Transform Global to Local (sub-pixel) Coordinates
    local_coords = list(range(0,nstacks))
     
    fit_x,fit_y,fit_z            = fit_params[:,2],fit_params[:,3],fit_params[:,4]
    global_x,global_y,global_z   = global_coords[:,0],global_coords[:,1],global_coords[:,2]
    mask_size_x, mask_size_y, mask_size_z = mask_size[0],mask_size[1],mask_size[2]
    
    centered_x = fit_x-mask_size_x/2
    centered_y = fit_y-mask_size_y/2
    
    if zsteps == 1:
        centered_z = 0
    else:
        centered_z = fit_z-mask_size_z/2
    
    local_x,local_y,local_z = global_x+centered_x, global_y+centered_y, global_z+centered_z
    local_coords_temp  = np.swapaxes(np.array([local_x,local_y,local_z],order='C'),0,1)
    local_coords = local_coords_temp
    
    return local_coords, fit_params, fit_errors

#%% Fit a 3D Gaussian to a ROI

def fit_3D_gauss(ROI_stack,mask_size_trace):
                 
    ##### Fit Data 
    # Get ROI-data from the stack
    data = np.ndarray.flatten(ROI_stack)
    size_x,size_y,size_z = mask_size_trace 
    x = np.arange(0,size_x)
    y = np.arange(0,size_y)
    z = np.arange(0,size_z)
    
    # Create xyz-coordinates 
    xx, yy, zz = np.meshgrid(x, y, z, sparse=False)
    xyz = [np.ndarray.flatten(xx),np.ndarray.flatten(yy), np.ndarray.flatten(zz)]      
    
    # Initial Fit Parameters
    in_A  = float(max(data))
    in_C  = float(min(data))
    in_x0 = max(x)/2
    in_y0 = max(y)/2
    in_z0 = max(z)/2
    in_wx = 1.5
    in_wy = 1.5
    in_wz = 3
    
    # p0 is the initial guess for the fit coefficients
    p0 = [in_A, in_C, in_x0, in_y0, in_z0, in_wx, in_wy, in_wz]
    
    try:
        coeff, var_matrix = curve_fit(gauss_3D(xyz,p0), xyz, data, p0=p0, absolute_sigma=True)
        perr = np.sqrt(np.diag(var_matrix))
        
    except RuntimeError:
        print('Error - curve_fit failed')
        coeff = np.empty(len(p0))                
        perr  = np.empty(len(p0))
        
    return coeff, perr

#%% Fit a 2D Gaussian to a ROI
#(Used when the measurement does not contain any z-stacks)
    
def fit_2D_gauss(ROI_stack,mask_size_trace):
                 
    ##### Fit Data 
    # Get ROI-data from the stack
    data = np.ndarray.flatten(ROI_stack)
    size_x,size_y,size_z = mask_size_trace 
    x = np.arange(0,size_x)
    y = np.arange(0,size_y)
    z = np.arange(0,size_z)
    
    # Create xyz-coordinates 
    xx, yy, zz = np.meshgrid(x, y, z, sparse=False)
    xyz = [np.ndarray.flatten(xx),np.ndarray.flatten(yy), np.ndarray.flatten(zz)]      
    
    # Initial Fit Parameters
    in_A  = float(max(data))
    in_C  = float(min(data))
    in_x0 = max(x)/2
    in_y0 = max(y)/2
    in_z0 = 0
    in_wx = 1.5
    in_wy = 1.5
    in_wz = 0
    
    # p0 is the initial guess for the fit coefficients
    p0 = [in_A, in_C, in_x0, in_y0, in_wx, in_wy]
    
    try:
        coeff, var_matrix = curve_fit(gauss_2D, xyz, data, p0=p0, absolute_sigma=True)
        perr = np.sqrt(np.diag(var_matrix))
        
    except RuntimeError:
        print('Error - curve_fit failed')
        coeff = np.empty(len(p0))                
        perr  = np.empty(len(p0))
        
    ## Workaround: inster z-values so .xlsx format stays in
    ## the same format as the 3D measurements     
    coeff = np.insert(coeff,4,in_z0)
    coeff = np.append(coeff,in_wz)
    
    perr  = np.insert(perr,4,0)
    perr  = np.append(perr,0)
    
    return coeff, perr

#%%
    
def fit_zprofile(file_path,global_coords,stack_nr,mask_size):
    # Get Measurement Parameters 
    file_log = file_io.read_logfile(file_path)
    zsteps = file_log[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])
    
    ### Get ROI from the stack using Global Coordinates and Fit 3D Gaussian

    # Make 3D Mask
    size_x,size_y,size_z=mask_size
    mask = create_3D_mask(size_x,size_y,size_z,10,10)
    
    # Define Bounds
    max_bounds=(65536,65536,size_x,size_y,size_z,size_x-1,size_y-1,size_z-1)
    #min_bounds=(0,0,0,0,0,0,0,0)
    #bounds = (min_bounds,max_bounds)
    
    ### Fit ROI to a 3D Gauss 
      
    # Read one stack from .bin file 
    if zsteps == 1:
        nstacks=1
        stack = file_io.read_bin_all(file_path)
    else:
        nframes = np.uint(str.split(file_log[7],'=')[1])
        nstacks = int(nframes/zsteps)
        stack = file_io.read_stack(file_path,stack_nr)
    
    # Allocate memory for fit parameters & errors for 1 stack 
    fit_params = np.empty([len(global_coords[:,0]),len(max_bounds)])
    fit_errors = np.empty([len(global_coords[:,0]),len(max_bounds)])
    mask_size_traces = np.empty([len(global_coords[:,0]),len(mask_size)])
    
    for trace_nr in range(0,len(global_coords)):
        
        # Get ROI
        peak_coords = global_coords[trace_nr,:]
        ROI_stack,mask_size_trace = func.get_ROI_from_stack(file_path,stack,peak_coords,mask)
        
        ##### Fit Data #####
        # Do not fit Z-coordinates when measurement is in 2D
        if zsteps ==1:
            coeff,perr = fit_2D_gauss(ROI_stack,mask_size_trace)
        else:
            coeff,perr = fit_3D_gauss(ROI_stack,mask_size_trace)
            
        fit_params[trace_nr,:]= coeff
        fit_errors[trace_nr,:]= perr
        mask_size_traces[trace_nr,:] = mask_size_trace
        

        
    ### Transform Global to Local (sub-pixel) Coordinates
    local_coords = list(range(0,nstacks))
     
    fit_x,fit_y,fit_z            = fit_params[:,2],fit_params[:,3],fit_params[:,4]
    global_x,global_y,global_z   = global_coords[:,0],global_coords[:,1],global_coords[:,2]
    mask_size_x, mask_size_y, mask_size_z = mask_size[0],mask_size[1],mask_size[2]
    
    centered_x = fit_x-mask_size_x/2
    centered_y = fit_y-mask_size_y/2
    
    if zsteps == 1:
        centered_z = 0
    else:
        centered_z = fit_z-mask_size_z/2
    
    local_x,local_y,local_z = global_x+centered_x, global_y+centered_y, global_z+centered_z
    local_coords_temp  = np.swapaxes(np.array([local_x,local_y,local_z],order='C'),0,1)
    local_coords = local_coords_temp
    
    return local_coords, fit_params, fit_errors