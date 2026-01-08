# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 13:40:36 2022

@author: grife

Finds the difference in point values between two files.
Difference for all arrays in file calculated unless list of array names given by 'array_name' parameter
A .vtk file with material id and differences between properties at idx

"""

import numpy as np
import vtk
from os.path import exists


def create_difference_plot(filename1, path_to_file_1, filename2, path_to_file_2, output_path, array_name=None):

    filename1 = filename1.split(".")[0]
    filename2 = filename2.split(".")[0]
    fullFilename1 = path_to_file_1 + filename1 + ".vtk"
    fullFilename2 = path_to_file_2 + filename2 + ".vtk"
    if exists(fullFilename1) and exists(fullFilename2):

        print("FILENAME 1: " + filename1)
        print("FILENAME 2: " + filename2)

        # Set up poly data reader for result set 1
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(fullFilename1)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        numArray = data.GetPointData().GetNumberOfArrays()

        # Set up poly data reader for result set 2
        reader2 = vtk.vtkUnstructuredGridReader()
        reader2.SetFileName(fullFilename2)
        reader2.ReadAllVectorsOn()
        reader2.ReadAllScalarsOn()
        reader2.Update()
        data2 = reader2.GetOutput()
        numArray2 = data2.GetPointData().GetNumberOfArrays()

        assert numArray == numArray2, numArray + " and  " + numArray2 + " are not equal"
        assert data.GetNumberOfPoints() == data2.GetNumberOfPoints(), data.GetNumberOfPoints() + " and  " + data2.GetNumberOfPoints() + " are not equal"

        # Create new unstrauctured grid
        output = vtk.vtkUnstructuredGridReader()
        output.SetFileName(fullFilename1)
        output.ReadAllVectorsOn()
        output.ReadAllScalarsOn()
        output.Update()
        output_data = data
        output_pointData = output_data.GetPointData()

        # Read in data and rename arrays
        for idx in range(0,numArray):
            arr_1 = data.GetPointData().GetArray(idx)
            arr_2 = data2.GetPointData().GetArray(idx)
            arr_full_name = arr_1.GetName()
            arr_full_name_2 = arr_2.GetName()
            ao = None
            ao1 = None
            ao2 = None
            if ("material" not in arr_full_name):
                if array_name is None or arr_full_name in array_name:

                    assert arr_full_name == arr_full_name_2, arr_full_name + " and  " + arr_full_name_2 + " are not the same"

                    # Rename array from file 2
                    prop_name1 = filename1 + "_" + arr_full_name
                    arr_1.SetName(prop_name1)

                    # Rename array from file 2
                    prop_name2 = filename2 + "_" + arr_full_name_2
                    arr_2.SetName(prop_name2)

                    # Calculate difference between two properties
                    difference = np.zeros(data2.GetNumberOfPoints())
                    arr_1_new = np.zeros(data2.GetNumberOfPoints())
                    arr_2_new = np.zeros(data2.GetNumberOfPoints())
                    for i in range(data2.GetNumberOfPoints()):
                        arr_1_val = arr_1.GetValue(i)
                        arr_2_val = arr_2.GetValue(i)
                        difference_in_vals = np.abs(arr_1_val - arr_2_val)
                        difference[i] = difference_in_vals
                        arr_1_new[i] = arr_1_val
                        arr_2_new[i] = arr_2_val

                    # Create new vector
                    ao = vtk.vtkFloatArray()
                    ao.SetName('Diff ' + arr_full_name)
                    ao.SetNumberOfValues(difference.size)
                    for x in range(difference.size):
                        ao.SetValue(x,difference[x])

                    # Update_vector
                    ao1 = vtk.vtkFloatArray()
                    ao1.SetName(prop_name1)
                    ao1.SetNumberOfValues(arr_1_new.size)
                    for x in range(arr_1_new.size):
                        ao1.SetValue(x,arr_1_new[x])

                    # Update_vector
                    ao2 = vtk.vtkFloatArray()
                    ao2.SetName(prop_name2)
                    ao2.SetNumberOfValues(arr_2_new.size)
                    for x in range(arr_2_new.size):
                        ao2.SetValue(x,arr_2_new[x])
            else:
                # Read in material_id array
                ao = data.GetPointData().GetArray(idx)
            if ao is not None:
                output_pointData.AddArray(ao)
                output_pointData.AddArray(ao1)
                output_pointData.AddArray(ao2)
                output_pointData.AddArray(ao1)
                output_pointData.AddArray(arr_1)
                output_pointData.AddArray(arr_2)
                print("Array " + ao.GetName() + " written")

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetInputData(output_data)
        file_output = filename1 + "_" + filename2 + ".vtk"
        output_full_filename = output_path + file_output
        writer.SetFileName(output_full_filename)
        writer.Update()
        writer.Write()

        print("COMPLETE: " + output_full_filename)
    else:
        print("Failed to open files, one or both of the following may not exists:")
        print(fullFilename1)
        print(fullFilename2)

if __name__ == "__main__":
    # path = "./output_atrophy/post_processing/OAS1_0005_MR1_atrophy_fine/"
    # file1 = "OAS1_0005_MR1_atrophy_19R_last_timestep"
    # file2 = "OAS1_0005_MR1_atrophy_1R_last_timestep"
    path = "./output_atrophy/post_processing/rampp_FA_9R/"
    file1 = "rampp_atrophy_9R_last_timestep"
    file2 = "rampp_atrophy_FA_last_timestep"
    out_path = path
    create_difference_plot(file1, path, file2, path, out_path)





