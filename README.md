# Readme for 2P_Specra

Overview of the Functions 

2P_spectra_main calls:

  - Get filepath
  - Get Global Coordinates
  - Fit Global Coordinates to fit_2DGaussian&(Lorentzian)^2 
  - Write all fit parameters and fit errors to .xlsx
  

functions_2PSpectra contains:

  - get_ROI_from_stack(filepath,stack,peak_coords,mask)
  - Get_Traces_3D(filepath,stack,threshold_factor,mask_size,max_num_peaks)
  - get_global_coords(filepath, threshold_factor, mask_size,max_num_peaks)
  - function_2Dgausslor_sq(xyl,*p)
  - fit_2Dgauss_lorsq(wavelength,ROI_stack,mask_size_trace)
  - get_local_and_spectrum(file_path,global_coords_corrected,mask_size,stack_nr)
  

file_io contains:

  - change_extension(path_file, extension):
  - get_path()
  - read_logfile(path_logfile)
  - read_bin(path_binfile,slice_nr)
  - read_stack(filepath,stack_nr)
  - read_bin_all(filepath)
  - write_xlsx_list(file_name,data_list,column_headers,column_indices)
  - write_xlsx_array(file_path,data_array,column_headers,column_indices)
  - write_pd2xlsx(file_path,pd_frame)
  

tools contains:
  - get_max(array):
  - get_min(array):
  - create_3D_mask(size_x,size_y,size_z,amp,sigma):
  - create_2D_mask(size_x,size_y,amp,sigma):
  - gauss_3D(xyz, p):
  - gauss_2D(xy, *p):
  - sq_lorentzian(l, *p):
  - get_zprofile(file_path,trace_coordinate,stack,mask_size):
  - fit_lorentzian_sq(x,data_z):
  - fit_spectra(file_path,global_coords,mask_size):
  - get_local_coordinates(file_path,global_coords,stack_nr,mask_size):
  - fit_3D_gauss(ROI_stack,mask_size_trace):
  - fit_2D_gauss(ROI_stack,mask_size_trace):
  - fit_zprofile(file_path,global_coords,stack_nr,mask_size):
  
