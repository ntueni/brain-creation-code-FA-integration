# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 14:41:55 2022

@author: grife

This creates a .csv file given the statistical information for each region within the files listed in 'fileNames'.
"""

import numpy as np
import vtk
import csv
from os.path import exists
from math import log10
import scipy.stats as sp
import matplotlib.pyplot as plt

path = "/path_to_post_processing_files"
fileNames = ["OAS1_0002_MR1_atrophy_19R_last_timestep","OAS1_0002_MR1_atrophy_9R_last_timestep","OAS1_0002_MR1_atrophy_4R_last_timestep","OAS1_0002_MR1_atrophy_2R_last_timestep", "OAS1_0002_MR1_atrophy_1R_last_timestep"]
out_path = path
filenameOUT = "OAS1_0002_MR1_atrophy_results.csv"
filenameDiff = "OAS1_0002_MR1_atrophy_diff.csv"
path_out = path + filenameOUT
path_diff = path + filenameDiff
dataMap = {}
measure = ['displacement', 'max_principal_stretch', 'med_principal_stretch', 'min_principal_stretch', 'max_shear', 'hyrdostatic_stress', 'von_mises']

for f in fileNames:
    written_material_id = False
    fullFilename1 = path + f + ".vtk"
    filenameData = {}
    if (exists(fullFilename1)):
        print("FILENAME: " + f)
        
        # Set up poly data reader for result set 1
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(fullFilename1)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        numArray = data.GetPointData().GetNumberOfArrays()
        
        # Read in data and rename arrays
        for idx in range(0, numArray):
            arr_1 = data.GetPointData().GetArray(idx)
            arr_full_name = arr_1.GetName()
            if "material" not in arr_full_name:
                name = arr_full_name.split("_")
                alt_arr_name = "_".join(name)
                
                py_arr_1 = [10000000]*data.GetNumberOfPoints()
                difference = np.zeros(data.GetNumberOfPoints())
                for i in range(data.GetNumberOfPoints()):
                    py_arr_1[i] = arr_1.GetValue(i)
                    
                # Save array data to map for file 1 and 2
                filenameData[alt_arr_name] = py_arr_1
                print("Array " + arr_full_name + " written")
                
            elif not written_material_id:
                written_material_id = True
                # Read in material_id array
                ao = data.GetPointData().GetArray(idx)
                mat_array = [10000000]*data.GetNumberOfPoints()
                for i in range(data.GetNumberOfPoints()):
                    mat_array[i] = ao.GetValue(i)
        
        # Add saved arrays/data to dataMap
        material_ids = mat_array
        region_data = {}
        for key, array in filenameData.items():
            region_values = {}
            for i in range(len(material_ids)):
                mat_id = material_ids[i]
                if region_values.__contains__(mat_id):
                    region_array = region_values[mat_id]
                else:
                    region_array = []
                region_array.append(array[i])
                region_values[mat_id] = region_array
            region_data[key] = region_values
        filenameData["regional_data"] = region_data
        dataMap[f] = filenameData

f_out = open(path_out, 'w', newline='')
writer = csv.writer(f_out)
header = ['Conditioned', 'Poisson ratio', 'ModelType', 'Quantity', 'Region', 'Location', 'Other Notes', 'n', 'mean', 'trimmed_mean_5', 'min',
          '5th percentile', '25th percentile', '50th percentile', '75th percentile', '95th percentile',
          'max', 'var', 'std_dev', 'trimmed_nobs', 'trimmed_mean', 'trimmed_min', 'trimmed_max', 'trimmed_var', 'trimmed_std_dev']
writer.writerow(header)
rows = []

for fileName, arrays in dataMap.items():
    [conditioning, poissons, modelType, *other] = fileName.split("_")
    location = other[0] if len(other) == 1 else 'L4'
    regional_data = arrays['regional_data']
    for mea, data_per_region in regional_data.items():
        if mea in measure:  # Check if mea is in measure
            quantity = mea
            for mat_id, mat_arrs in data_per_region.items():
                nobs = len(mat_arrs)
                mean = sp.tmean(mat_arrs)
                trimmed_mean_5 = sp.trim_mean(mat_arrs, 0.05)
                min_value = min(mat_arrs)
                perc_5th = sp.scoreatpercentile(mat_arrs, 5)
                perc_25th = sp.scoreatpercentile(mat_arrs, 25)
                perc_50th = sp.scoreatpercentile(mat_arrs, 50)
                perc_75th = sp.scoreatpercentile(mat_arrs, 75)
                perc_95th = sp.scoreatpercentile(mat_arrs, 95)
                max_value = max(mat_arrs)
                variance = sp.tvar(mat_arrs)
                std_dev = sp.tstd(mat_arrs)
                
                iqr = sp.iqr(mat_arrs)
                lower_limit = perc_25th - (1.5 * iqr)
                upper_limit = perc_75th + (1.5 * iqr)
                lower_idx = 0
                upper_idx = len(mat_arrs)
                mat_arrs.sort()
                for idx, i in enumerate(mat_arrs):
                    if i > lower_limit:
                        lower_idx = idx
                        break
                for idx in range(len(mat_arrs) - 1, -1, -1):
                    i = mat_arrs[idx]
                    if i < upper_limit:
                        upper_idx = idx
                        break
                
                trimmed_mat = mat_arrs[lower_idx:upper_idx + 1]
                trimmed_nobs = len(trimmed_mat)
                trimmed_mean = sp.tmean(trimmed_mat)
                trimmed_min_value = min(trimmed_mat)
                trimmed_max_value = max(trimmed_mat)
                trimmed_variance = sp.tvar(trimmed_mat)
                trimmed_std_dev = sp.tstd(trimmed_mat)
                row = [conditioning, poissons, modelType, quantity, mat_id, location, "_".join(other),
                       nobs, mean, trimmed_mean_5, min_value, perc_5th, perc_25th, perc_50th, perc_75th, perc_95th, max_value, variance, std_dev,
                       trimmed_nobs, trimmed_mean, trimmed_min_value, trimmed_max_value, trimmed_variance, trimmed_std_dev]

                rows.append(row)

writer.writerows(rows)

f_out.close()

