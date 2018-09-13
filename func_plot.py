# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 15:28:13 2018

@author: Redmar
"""

#%% Plot-related Functions for 3D-Tracking data

import numpy as np
import functions as func
import file_io as file_io
import matplotlib.pyplot as plt


#%% Plot all found cooridnates in a 3D plot (X,Y,Z)
def plot_3D(coords,fig_num,filepath):
    
    logfile = file_io.read_logfile(filepath)

    ypix = np.uint(str.split(logfile[8],'=')[1])
    xpix = np.uint(str.split(logfile[9],'=')[1])
    zsteps = logfile[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0]) 
    
    num_stacks = len(coords)
    for stack_nr in range(0,num_stacks):
        
        x_loc = coords[stack_nr][:,0]
        y_loc = coords[stack_nr][:,1]
        z_loc = coords[stack_nr][:,2]
    
        fig = plt.figure(fig_num) 
        ax = fig.gca(projection='3d')
        ax.scatter(x_loc,y_loc,z_loc,zdir='z')
    
    plt.axis([0,xpix,0,ypix,0,zsteps])
    plt.xlabel('X-Coordinate [pix]')
    plt.ylabel('Y-Coorindate [pix]')
    return logfile

#%% Plot Histogram of Fit Errors 
def plot_hist_errs(errors,fig_num,filepath):
    
    logfile = file_io.read_logfile(filepath)
    
    zsteps = str.split(logfile[19],'=')[1]
    zsteps = zsteps.replace(',','.',1)
    zsteps = np.double(zsteps)
    
    pix_size = str.split(logfile[10],'=')[1]
    pix_size = pix_size.replace(',','.',1)
    pix_size = np.double(pix_size)
    pix_size = (pix_size*260)/pix_size
    
    fit_errors = errors
    
    plt.figure(fig_num)
    plt.title('Fit Errors')
    
    fit_err_x = np.ndarray.flatten(fit_errors[0][:,2])*pix_size
    fit_err_y = np.ndarray.flatten(fit_errors[0][:,3])*pix_size
    fit_err_z = np.ndarray.flatten(fit_errors[0][:,4])*zsteps
    
    plt.subplot(3,1,1)
    plt.title('X-Coordinate',fontsize=16)
    plt.hist(fit_err_x,50,[0,10])
    plt.tick_params(labelsize=16)
    
    plt.subplot(3,1,2)
    plt.title('Y-Coordinate',fontsize=16)
    plt.hist(fit_err_y,50,[0,10])
    plt.tick_params(labelsize=16)
    
    plt.subplot(3,1,3)
    plt.title('Z-Coordinate',fontsize=16)
    plt.hist(fit_err_z,50,[0,10])
    plt.tight_layout()
    
    plt.xlabel('Fit Error [nm]',fontsize=16)
    plt.tick_params(labelsize=16)
    
    
#%% PLot a single trace  
def plot_single_trace(stacks_sorted,trace_num,fig_num):

    for i in range(0,len(stacks_sorted)):
        try:
            x,y,z = stacks_sorted[i][trace_num,0],stacks_sorted[i][trace_num,1],stacks_sorted[i][trace_num,2]
            fig = plt.figure(fig_num) 
            ax = fig.gca(projection='3d')
            ax.scatter(x,y,z,zdir='z')
            
            
        except IndexError:
            print(':-(')

#%% Plot Raw trace data with fitted Gaussian (3D)        
def plot_data_w_fit(filepath,global_coords,peak,fit_params): 
    mask = func.create_3D_mask(11,11,23,10,10)
    stack = file_io.get_stack(filepath,0)
    peak=81
    ROI_stack,mask_size = func.get_ROI_from_stack(filepath,stack,global_coords[0][peak],mask)
    
    data = np.ndarray.flatten(ROI_stack)
    
    size_x,size_y,size_z = mask_size
    x = np.arange(0,size_x)
    y = np.arange(0,size_y)
    z = np.arange(0,size_z)
    
    # Create xyz-coordinates 
    xx, yy, zz = np.meshgrid(x, y, z, sparse=False)
    xyz = [np.ndarray.flatten(xx),np.ndarray.flatten(yy), np.ndarray.flatten(zz)]
    params = fit_params[0][peak]
    fit = func.gauss_3D(xyz,*params)
    
    data = data-np.median(data)
    fit  = fit -np.median(fit)
    
    plt.figure(5)
    line_data = plt.plot(data,'.r' )
    plt.setp(line_data, 'linewidth',5)
    line_fit = plt.plot(fit,'-k')
    plt.setp(line_fit, 'linewidth',1)
    plt.axis([0,len(fit),-5, max(data)+10])
    plt.xlabel('Pixel Index [#]',fontsize=16)
    plt.ylabel('Intensity [a.u.]',fontsize=16)
    plt.tick_params(labelsize=16)