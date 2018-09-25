# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 14:57:53 2018

@author: vlieg
"""
import numpy as np
import tkinter as tk
import file_io as file_io
from scipy.optimize import curve_fit
from tkinter import filedialog
import pandas as pd 
import tools as tools 

#%% Get ROI from stack 

def get_ROI_from_stack(filepath,stack,peak_coords,mask):
    
    # Get file parameters     
    logfile = file_io.read_logfile(filepath)

    ypix = np.uint(str.split(logfile[8],'=')[1])
    xpix = np.uint(str.split(logfile[9],'=')[1])
    zsteps = logfile[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])   
    size_x,size_y,size_z=mask.shape
    
    x_range = np.int16([peak_coords[0]-np.floor(size_x/2),peak_coords[0]+np.floor(size_x/2)])
    y_range = np.int16([peak_coords[1]-np.floor(size_y/2),peak_coords[1]+np.floor(size_y/2)])
    z_range = np.int16([peak_coords[2]-np.floor(size_z/2),peak_coords[2]+np.floor(size_z/2)])
    
    mask_bounds = np.int16([np.floor(size_x/2),np.floor(size_y/2), np.floor(size_z/2)])    
    mask_temp = mask
    
    # X-coordinates too small:
    if  peak_coords[0] < mask_bounds[0]:
        x_range = np.int16([0,peak_coords[0]+np.floor(size_x/2)])
        mask_temp = mask_temp[mask_bounds[0]-peak_coords[0]::,:,:]
     
    # X-coordinates too large:
    if peak_coords[0] + mask_bounds[0] >= xpix:
        x_range = np.int16([peak_coords[0]-np.floor(size_x/2),xpix])
        ind_cut = ((x_range[0])+size_x)-xpix
        mask_temp = mask_temp[0:(size_x-ind_cut),:,:]
    
    # Y-coordinates too small:
    if peak_coords[1] < mask_bounds[1]:
        y_range = np.int16([0,peak_coords[1]+np.floor(size_y/2)])
        mask_temp = mask_temp[:,mask_bounds[1]-peak_coords[1]::,:]
        
    # Y-coordinates too large:             
    if peak_coords[1] + mask_bounds[1] >= ypix:
        y_range = np.int16([peak_coords[1]-np.floor(size_y/2),ypix])
        ind_cut = ((y_range[0])+size_y)-ypix
        mask_temp = mask_temp[:,0:(size_y-ind_cut),:]
    
    # Z-coordinates too small:          
    if peak_coords[2] < mask_bounds[2] and zsteps != 1:
        z_range = np.int16([0,peak_coords[2]+np.floor(size_z/2)])
        mask_temp = mask_temp[:,:,mask_bounds[2]-peak_coords[2]::]
        
    # Z-coordinates too large:
    if peak_coords[2] + mask_bounds[2] >= zsteps:
        z_range = np.int16([peak_coords[2]-np.floor(size_z/2),zsteps])
        ind_cut = ((z_range[0])+size_z)-zsteps
        mask_temp = mask_temp[:,:,0:(size_z-ind_cut)]
    
    if zsteps != 1:    
        ROI_stack = stack[x_range[0]:x_range[1]+1,y_range[0]:y_range[1]+1,z_range[0]:z_range[1]+1]
        mask_size = np.shape(mask_temp)
    
    else:
        ROI_stack = stack[x_range[0]:x_range[1]+1,y_range[0]:y_range[1]+1,:]
        mask_size = np.shape(mask_temp)
        
    return ROI_stack, mask_size
    
    
#%% Get Trace Coordinates (X,Y,Z) from a stack
    
def Get_Traces_3D(filepath,stack,threshold_factor,mask_size,max_num_peaks):
        
    # Get file parameters     
    logfile = file_io.read_logfile(filepath)

    ypix = np.uint(str.split(logfile[8],'=')[1])
    xpix = np.uint(str.split(logfile[9],'=')[1])
    zsteps = logfile[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])    
        
    # Make 3D Mask
    size_x, size_y, size_z = mask_size 
    mask_bounds = np.int16([np.floor(size_x/2),np.floor(size_y/2), np.floor(size_z/2)])
    mask = tools.create_3D_mask(size_x, size_y, size_z,10,10)
    
    # Find maximum intensity indeces
    
    threshold = np.median(stack)*threshold_factor
    max_intensity=np.max(stack)
    masked_stack = stack
    peak_coordinates = np.array([], dtype=np.uint16)
    num_peaks = 0
    
    if not max_num_peaks or max_num_peaks is 0:
        max_num_peaks = np.inf
    
    while max_intensity > threshold and num_peaks <= max_num_peaks:
    
        peak_coords = np.unravel_index(np.argmax(masked_stack), np.shape(masked_stack))
        max_intensity = masked_stack[peak_coords]
        peak_coords = np.int16(peak_coords)    
        
        x_range = np.int16([peak_coords[0]-np.floor(size_x/2),peak_coords[0]+np.floor(size_x/2)])
        y_range = np.int16([peak_coords[1]-np.floor(size_y/2),peak_coords[1]+np.floor(size_y/2)])
        z_range = np.int16([peak_coords[2]-np.floor(size_z/2),peak_coords[2]+np.floor(size_z/2)])
        
        mask_temp = mask
    
        # X-coordinates too small:
        if  peak_coords[0] < mask_bounds[0]:
            x_range = np.int16([0,peak_coords[0]+np.floor(size_x/2)])
            mask_temp = mask_temp[mask_bounds[0]-peak_coords[0]::,:,:]
         
        # X-coordinates too large:
        if peak_coords[0] + mask_bounds[0] >= xpix:
            x_range = np.int16([peak_coords[0]-np.floor(size_x/2),xpix])
            ind_cut = ((x_range[0])+size_x)-xpix
            mask_temp = mask_temp[0:(size_x-ind_cut),:,:]
    
        # Y-coordinates too small:
        if peak_coords[1] < mask_bounds[1]:
            y_range = np.int16([0,peak_coords[1]+np.floor(size_y/2)])
            mask_temp = mask_temp[:,mask_bounds[1]-peak_coords[1]::,:]
            
        # Y-coordinates too large:             
        if peak_coords[1] + mask_bounds[1] >= ypix:
            y_range = np.int16([peak_coords[1]-np.floor(size_y/2),ypix])
            ind_cut = ((y_range[0])+size_y)-ypix
            mask_temp = mask_temp[:,0:(size_y-ind_cut),:]
            
        # Z-coordinates too small:          
        if peak_coords[2] < mask_bounds[2]:
            z_range = np.int16([0,peak_coords[2]+np.floor(size_z/2)])
            mask_temp = mask_temp[:,:,mask_bounds[2]-peak_coords[2]::]
            
        # Z-coordinates too large:
        if peak_coords[2] + mask_bounds[2] >= zsteps:
            z_range = np.int16([peak_coords[2]-np.floor(size_z/2),zsteps])
            ind_cut = ((z_range[0])+size_z)-zsteps
            mask_temp = mask_temp[:,:,0:(size_z-ind_cut)]
               
        ROI_stack = stack[x_range[0]:x_range[1]+1,y_range[0]:y_range[1]+1,z_range[0]:z_range[1]+1]
        ROI_stack_masked = ROI_stack*mask_temp
        masked_stack[x_range[0]:x_range[1]+1,y_range[0]:y_range[1]+1,z_range[0]:z_range[1]+1]=ROI_stack_masked
        peak_coordinates = np.append(peak_coordinates,peak_coords)
        
        num_peaks = num_peaks+1
    
    peak_coordinates = np.reshape(peak_coordinates,[3,int(len(peak_coordinates)/3)],1)
    peak_coordinates = np.transpose(peak_coordinates)
    
    return peak_coordinates, num_peaks
    
    
#%% Get the global coordinates of traces from a 3D stack of images 
    
def get_global_coords(filepath, threshold_factor, mask_size,max_num_peaks):
    # Get file parameters 
    logfile = file_io.read_logfile(filepath)
    
    ypix   = np.uint(str.split(logfile[8],'=')[1])
    xpix   = np.uint(str.split(logfile[9],'=')[1])
    zsteps = logfile[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])
    nframes= np.uint(str.split(logfile[7],'=')[1])
    nstacks= int(nframes/zsteps)
    
    if zsteps == 1:
        nstacks=1
             
    else:
        nframes = np.uint(str.split(logfile[7],'=')[1])
        nstacks = int(nframes/zsteps)
       
    # Get Peak coordinates from all stacks 
    stack=np.zeros([xpix,ypix,zsteps], dtype = np.uint16)
    num_traces_pstack=np.empty(nstacks, dtype = np.uint16)
    peak_coordinates = list(range(0,nstacks))
    
    for stack_nr in range(0,nstacks):
        # Read stack from file     
        for slice_nr in range(stack_nr*zsteps,stack_nr*zsteps + zsteps):
            stack[:,:,slice_nr-stack_nr*zsteps]=file_io.read_bin(filepath,slice_nr)
            
        # Get Peak coordinates from stack
        [peak_coordinates_stack, num_trace] = Get_Traces_3D(filepath,stack,threshold_factor,mask_size,max_num_peaks)
        peak_coordinates[stack_nr] = peak_coordinates_stack
        num_traces_pstack[stack_nr]=num_trace
    
    ## Remove traces with wrong coordinates
    if len(peak_coordinates) == 1:
        peak_coordinates=peak_coordinates[0]
        
    global_coords_corrected = np.empty([0,3])
    trace_nr=0
    for trace_coords in peak_coordinates:
        if trace_coords[0]<5 or trace_coords[1]<5:
            trace_nr=trace_nr+1
        else:
            trace_coords = np.expand_dims(trace_coords,1)
            trace_coords = np.swapaxes(trace_coords,0,1)
            global_coords_corrected=np.append(global_coords_corrected,trace_coords,0)
            trace_nr=trace_nr+1
            
    return global_coords_corrected, num_traces_pstack
    


#%% Fit Function for x,y,lambda 
    
def function_2Dgausslor_sq(xyl,*p):
    xx,yy,ll = xyl[0],xyl[1],xyl[2]
    A,C,x0,y0,wx,wy,l0,gamma = p
    
    gauss_x = np.square(xx-x0)/(2*np.square(wx))
    gauss_y = np.square(yy-y0)/(2*np.square(wy))
    lorentz = 1/(1+4*(np.sqrt(2)-1)*np.square((ll-l0)/gamma))
    
    return (A*np.exp(-1*(gauss_x+gauss_y))*lorentz)+C


#%%
    
def fit_2Dgauss_lorsq(wavelength,ROI_stack,mask_size_trace):
                 
        ##### Fit Data #####
    # Get ROI-data from the stack
    data = np.ndarray.flatten(ROI_stack,order='F')

    xx,yy,ll  = np.meshgrid(np.arange(mask_size_trace[0]),np.arange(0,mask_size_trace[1]),wavelength)
    xyl     = [np.ndarray.flatten(yy,order='F'),np.ndarray.flatten(xx,order='F'),np.ndarray.flatten(ll,order='F')]
    
    # Initial Fit Parameters
    in_A  = float(max(data))
    in_C  = float(min(data))
    in_x0 = mask_size_trace[0]/2
    in_y0 = mask_size_trace[1]/2
    in_wx = 1.5
    in_wy = 1.5
    in_l0 = xyl[2][tools.get_max(data)[1]]
    in_gamma = 35
    
    # p0 is the initial guess for the fit coefficients
    p0 = [in_A, in_C, in_x0, in_y0, in_wx, in_wy, in_l0, in_gamma]
    
    try:
        coeff, var_matrix = curve_fit(function_2Dgausslor_sq, xyl, data, p0=p0, absolute_sigma=True)
        perr = np.sqrt(np.diag(var_matrix))
        
        #plt.figure(i)
        #plt.plot(np.arange(0,len(data)),data,np.arange(0,len(data)),func.function_2Dgausslor_sq(xyl,*coeff))
        #plt.title(['Trace# =',str(i)])
    except RuntimeError:
        print('Error - curve_fit failed')
        coeff = np.empty(len(p0))                
        perr  = np.empty(len(p0))

    
    return coeff, perr


#%% Fit VOI to a 2D Gauss & Lorentzian
    
def get_local_and_spectrum(file_path,global_coords_corrected,mask_size,stack_nr):
    i=0
    file_path_dat = file_io.change_extension(file_path,'.dat')
    data_dat = pd.read_csv(file_path_dat,delimiter='\t',decimal =',')
    wavelength = np.array(data_dat['Spectrum peak (nm)'])
    
    # Get Measurement Parameters 
    file_log = file_io.read_logfile(file_path)
    zsteps = file_log[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])
    
    ### Get ROI from the stack using Global Coordinates and Fit 3D Gaussian

    # Make 3D Mask
    size_x,size_y,size_z=mask_size
    
    # Define Bounds
    max_bounds=(65536,65536,size_x,size_y,size_x-1,size_y-1,1200,300)
    #min_bounds=(0,0,0,0,0,0,0,0)
    #bounds = (min_bounds,max_bounds)
    
    ### Fit ROI to a 3D Gauss 
      
    # Read one stack from .bin file 
    if zsteps == 1:
        nstacks = 1
        stack   = file_io.read_bin_all(file_path)
        mask = tools.create_3D_mask(size_x,size_y,len(wavelength),10,10)
        
    else:
        nframes = np.uint(str.split(file_log[7],'=')[1])
        nstacks = int(nframes/zsteps)
        stack = file_io.read_stack(file_path,stack_nr)
        mask = tools.create_3D_mask(size_x,size_y,size_z,10,10)
    
    # Allocate memory for fit parameters & errors for 1 stack 
    fit_params = np.empty([len(global_coords_corrected[:,0]),len(max_bounds)])
    fit_errors = np.empty([len(global_coords_corrected[:,0]),len(max_bounds)])
    mask_size_traces = np.empty([len(global_coords_corrected[:,0]),len(mask_size)])
    
    
    for trace_nr in range(0,len(global_coords_corrected)):

        # Get ROI
        peak_coords = global_coords_corrected[trace_nr,:]
        ROI_stack,mask_size_trace = get_ROI_from_stack(file_path,stack,peak_coords,mask)
        
        ##### Fit Data #####
        # Do not fit Z-coordinates when measurement is in 2D
        
        if zsteps ==1:
            coeff,perr = fit_2Dgauss_lorsq(wavelength,ROI_stack,mask_size_trace)        
            i=i+1
            print(i)   
            
        else:
            coeff,perr = tools.fit_3D_gauss(ROI_stack,mask_size_trace)
            
            
        fit_params[trace_nr,:]= coeff
        fit_errors[trace_nr,:]= perr
        mask_size_traces[trace_nr,:] = mask_size_trace
        

    traces_parameters = np.array(np.concatenate((global_coords_corrected,fit_params,fit_errors),1))


    column_headers= ['x_global (pix)',
                 'y_global (pix)',
                 'z_global (pix)',
                 'amplitude (a.u.)',
                 'offset (a.u.)',
                 'x_local (pix)',
                 'y_local (pix)',
                 'width_x (pix)',
                 'width_y (pix)',
                 'spectrum peak (nm)',
                 'spectrum width (nm)',
                 'err_amplitude (a.u.)',
                 'err_offset (a.u.)',
                 'err x_local (pix)',
                 'err_y_local (pix)',
                 'err_width_x (pix)', 
                 'err_width_y (pix)', 
                 'err_spectrum peak (nm)',
                 'err_peak width (nm)']

    ## Convert Fit Parameters to DataFrame
    pd_fitparams = pd.DataFrame(data=traces_parameters,columns=column_headers)
        

    ### Transform Global to Local (sub-pixel) Coordinates
    local_coords = list(range(0,nstacks))
     
    fit_x,fit_y,fit_z            = fit_params[:,2],fit_params[:,3],fit_params[:,4]
    global_x,global_y,global_z   = global_coords_corrected[:,0],global_coords_corrected[:,1],global_coords_corrected[:,2]
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
    
    pd_fitparams['x_local (pix)'] = local_coords[:,0]
    pd_fitparams['y_local (pix)'] = local_coords[:,1]
    
    return pd_fitparams