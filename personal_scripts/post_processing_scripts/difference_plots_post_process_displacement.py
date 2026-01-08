# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 13:40:36 2022
Modified to compute displacement magnitude differences

@author: grife

Finds the difference in displacement magnitudes between two files.
Computes displacement magnitude from 3 components (x, y, z) and calculates differences.
A .vtk file with material id and differences between displacement magnitudes.

"""

import numpy as np
import vtk
from os.path import exists


def create_displacement_magnitude_difference_plot(filename1, path_to_file_1, filename2, path_to_file_2, output_path, 
                                                displacement_array_name="displacement"):
    """
    Creates a difference plot between displacement magnitudes from two VTK files.
    
    Parameters:
    -----------
    filename1, filename2 : str
        Names of the VTK files (with or without .vtk extension)
    path_to_file_1, path_to_file_2 : str
        Paths to the respective files
    output_path : str
        Path where the output file will be saved
    displacement_array_name : str
        Name of the displacement vector field array in the VTK files
    """

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

        # Set up poly data reader for result set 2
        reader2 = vtk.vtkUnstructuredGridReader()
        reader2.SetFileName(fullFilename2)
        reader2.ReadAllVectorsOn()
        reader2.ReadAllScalarsOn()
        reader2.Update()
        data2 = reader2.GetOutput()

        # Verify both datasets have same number of points
        assert data.GetNumberOfPoints() == data2.GetNumberOfPoints(), \
            f"{data.GetNumberOfPoints()} and {data2.GetNumberOfPoints()} are not equal"

        # Get displacement arrays from both files
        disp_array_1 = data.GetPointData().GetArray(displacement_array_name)
        disp_array_2 = data2.GetPointData().GetArray(displacement_array_name)
        
        if disp_array_1 is None:
            print(f"Warning: Displacement array '{displacement_array_name}' not found in file 1")
            return
        if disp_array_2 is None:
            print(f"Warning: Displacement array '{displacement_array_name}' not found in file 2")
            return

        num_points = data.GetNumberOfPoints()
        
        # Calculate displacement magnitudes for both datasets
        displacement_mag_1 = np.zeros(num_points)
        displacement_mag_2 = np.zeros(num_points)
        displacement_mag_diff = np.zeros(num_points)
        
        print("Computing displacement magnitudes...")
        for i in range(num_points):
            # Get displacement components for point i from both files
            disp_1 = disp_array_1.GetTuple3(i)  # Returns (x, y, z) tuple
            disp_2 = disp_array_2.GetTuple3(i)  # Returns (x, y, z) tuple
            
            # Calculate magnitudes: |disp| = sqrt(x² + y² + z²)
            mag_1 = np.sqrt(disp_1[0]**2 + disp_1[1]**2 + disp_1[2]**2)
            mag_2 = np.sqrt(disp_2[0]**2 + disp_2[1]**2 + disp_2[2]**2)
            
            displacement_mag_1[i] = mag_1
            displacement_mag_2[i] = mag_2
            displacement_mag_diff[i] = np.abs(mag_1 - mag_2)

        # Create output data structure based on first file
        output_data = vtk.vtkUnstructuredGrid()
        output_data.DeepCopy(data)
        output_pointData = output_data.GetPointData()
        
        # Clear existing arrays to avoid duplication
        output_pointData.Initialize()

        # Add material_ids array if it exists (preserve as-is, no difference calculation)
        material_array = data.GetPointData().GetArray("material_ids")
        if material_array is not None:
            output_pointData.AddArray(material_array)
            print("Material_ids array preserved")
        else:
            # Try alternative name "material" as fallback
            material_array = data.GetPointData().GetArray("material")
            if material_array is not None:
                output_pointData.AddArray(material_array)
                print("Material array preserved (found as 'material')")

        # Create VTK arrays for displacement magnitudes and differences
        # Displacement magnitude from file 1
        mag_array_1 = vtk.vtkFloatArray()
        mag_array_1.SetName(f"{filename1}_displacement_magnitude")
        mag_array_1.SetNumberOfValues(num_points)
        for i in range(num_points):
            mag_array_1.SetValue(i, displacement_mag_1[i])

        # Displacement magnitude from file 2
        mag_array_2 = vtk.vtkFloatArray()
        mag_array_2.SetName(f"{filename2}_displacement_magnitude")
        mag_array_2.SetNumberOfValues(num_points)
        for i in range(num_points):
            mag_array_2.SetValue(i, displacement_mag_2[i])

        # Difference in displacement magnitudes
        diff_array = vtk.vtkFloatArray()
        diff_array.SetName("displacement_magnitude_difference")
        diff_array.SetNumberOfValues(num_points)
        for i in range(num_points):
            diff_array.SetValue(i, displacement_mag_diff[i])

        # Add arrays to output - focusing on the difference
        output_pointData.AddArray(mag_array_1)
        output_pointData.AddArray(mag_array_2)
        output_pointData.AddArray(diff_array)

        # Set the difference array as the active scalar for visualization
        output_pointData.SetActiveScalars("displacement_magnitude_difference")
        
        print(f"\nOutput arrays created:")
        print(f"  - {filename1}_displacement_magnitude")
        print(f"  - {filename2}_displacement_magnitude") 
        print(f"  - displacement_magnitude_difference (ACTIVE SCALAR)")
        print(f"  - Active scalar set to displacement_magnitude_difference for visualization")

        # Write output file
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetInputData(output_data)
        file_output = f"{filename1}_{filename2}_displacement_magnitude_diff.vtk"
        output_full_filename = output_path + file_output
        writer.SetFileName(output_full_filename)
        writer.Update()
        writer.Write()

        # Calculate percentiles for outlier removal
        p5 = np.percentile(displacement_mag_diff, 5)
        p95 = np.percentile(displacement_mag_diff, 95)
        
        # Filter out outliers (keep values between 5th and 95th percentiles)
        filtered_diff = displacement_mag_diff[(displacement_mag_diff >= p5) & (displacement_mag_diff <= p95)]
        
        # Print statistics
        print(f"\nDisplacement magnitude statistics:")
        print(f"File 1 ({filename1}):")
        print(f"  Min: {np.min(displacement_mag_1):.6f}")
        print(f"  Max: {np.max(displacement_mag_1):.6f}")
        print(f"  Mean: {np.mean(displacement_mag_1):.6f}")
        
        print(f"File 2 ({filename2}):")
        print(f"  Min: {np.min(displacement_mag_2):.6f}")
        print(f"  Max: {np.max(displacement_mag_2):.6f}")
        print(f"  Mean: {np.mean(displacement_mag_2):.6f}")
        
        print(f"Magnitude Difference (all values):")
        print(f"  Min: {np.min(displacement_mag_diff):.6f}")
        print(f"  Max: {np.max(displacement_mag_diff):.6f}")
        print(f"  Mean: {np.mean(displacement_mag_diff):.6f}")
        
        print(f"Magnitude Difference (outliers removed - 5th to 95th percentiles):")
        print(f"  5th percentile: {p5:.6f}")
        print(f"  95th percentile: {p95:.6f}")
        print(f"  Min (after outlier removal): {np.min(filtered_diff):.6f}")
        print(f"  Max (after outlier removal): {np.max(filtered_diff):.6f}")
        print(f"  Mean (after outlier removal): {np.mean(filtered_diff):.6f}")
        print(f"  Number of points used: {len(filtered_diff)}/{len(displacement_mag_diff)} ({100*len(filtered_diff)/len(displacement_mag_diff):.1f}%)")

        print(f"\nCOMPLETE: {output_full_filename}")
        
    else:
        print("Failed to open files, one or both of the following may not exist:")
        print(fullFilename1)
        print(fullFilename2)


# def create_displacement_magnitude_difference_plot(filename1, path_to_file_1, filename2, path_to_file_2, output_path, 
#                                                 displacement_array_name="displacement"):
#     """
#     Creates a difference plot between displacement magnitudes from two VTK files.
    
#     Parameters:
#     -----------
#     filename1, filename2 : str
#         Names of the VTK files (with or without .vtk extension)
#     path_to_file_1, path_to_file_2 : str
#         Paths to the respective files
#     output_path : str
#         Path where the output file will be saved
#     displacement_array_name : str
#         Name of the displacement vector field array in the VTK files
#     """

#     filename1 = filename1.split(".")[0]
#     filename2 = filename2.split(".")[0]
#     fullFilename1 = path_to_file_1 + filename1 + ".vtk"
#     fullFilename2 = path_to_file_2 + filename2 + ".vtk"
    
#     if exists(fullFilename1) and exists(fullFilename2):

#         print("FILENAME 1: " + filename1)
#         print("FILENAME 2: " + filename2)

#         # Set up poly data reader for result set 1
#         reader = vtk.vtkUnstructuredGridReader()
#         reader.SetFileName(fullFilename1)
#         reader.ReadAllVectorsOn()
#         reader.ReadAllScalarsOn()
#         reader.Update()
#         data = reader.GetOutput()

#         # Set up poly data reader for result set 2
#         reader2 = vtk.vtkUnstructuredGridReader()
#         reader2.SetFileName(fullFilename2)
#         reader2.ReadAllVectorsOn()
#         reader2.ReadAllScalarsOn()
#         reader2.Update()
#         data2 = reader2.GetOutput()

#         # Verify both datasets have same number of points
#         assert data.GetNumberOfPoints() == data2.GetNumberOfPoints(), \
#             f"{data.GetNumberOfPoints()} and {data2.GetNumberOfPoints()} are not equal"

#         # Get displacement arrays from both files
#         disp_array_1 = data.GetPointData().GetArray(displacement_array_name)
#         disp_array_2 = data2.GetPointData().GetArray(displacement_array_name)
        
#         if disp_array_1 is None:
#             print(f"Warning: Displacement array '{displacement_array_name}' not found in file 1")
#             return
#         if disp_array_2 is None:
#             print(f"Warning: Displacement array '{displacement_array_name}' not found in file 2")
#             return

#         num_points = data.GetNumberOfPoints()
        
#         # Calculate displacement magnitudes for both datasets
#         displacement_mag_1 = np.zeros(num_points)
#         displacement_mag_2 = np.zeros(num_points)
#         displacement_mag_diff = np.zeros(num_points)
#         signed_diff_array = np.zeros(num_points)
#         relative_diff_array = np.zeros(num_points)
        
#         print("Computing displacement magnitudes and differences...")
#         for i in range(num_points):
#             # Get displacement components for point i from both files
#             disp_1 = disp_array_1.GetTuple3(i)  # Returns (x, y, z) tuple
#             disp_2 = disp_array_2.GetTuple3(i)  # Returns (x, y, z) tuple
            
#             # Calculate magnitudes: |disp| = sqrt(x² + y² + z²)
#             mag_1 = np.sqrt(disp_1[0]**2 + disp_1[1]**2 + disp_1[2]**2)
#             mag_2 = np.sqrt(disp_2[0]**2 + disp_2[1]**2 + disp_2[2]**2)
            
#             displacement_mag_1[i] = mag_1
#             displacement_mag_2[i] = mag_2
            
#             # Compute signed and absolute differences
#             signed_diff = mag_1 - mag_2
#             displacement_mag_diff[i] = np.abs(signed_diff)
#             signed_diff_array[i] = signed_diff
            
#             # Relative difference (useful for comparing different scales)
#             avg_mag = (mag_1 + mag_2) / 2.0
#             if avg_mag > 1e-12:  # Avoid division by very small numbers
#                 relative_diff_array[i] = signed_diff / avg_mag
#             else:
#                 relative_diff_array[i] = 0.0

#         # Create output data structure based on first file
#         output_data = vtk.vtkUnstructuredGrid()
#         output_data.DeepCopy(data)
#         output_pointData = output_data.GetPointData()
        
#         # Clear existing arrays to avoid duplication
#         output_pointData.Initialize()

#         # Add material_ids array if it exists (preserve as-is, no difference calculation)
#         material_array = data.GetPointData().GetArray("material_ids")
#         if material_array is not None:
#             output_pointData.AddArray(material_array)
#             print("Material_ids array preserved")
#         else:
#             # Try alternative name "material" as fallback
#             material_array = data.GetPointData().GetArray("material")
#             if material_array is not None:
#                 output_pointData.AddArray(material_array)
#                 print("Material array preserved (found as 'material')")

#         # Create VTK arrays for displacement magnitudes and differences
#         # Displacement magnitude from file 1
#         mag_array_1 = vtk.vtkFloatArray()
#         mag_array_1.SetName(f"{filename1}_displacement_magnitude")
#         mag_array_1.SetNumberOfValues(num_points)
#         for i in range(num_points):
#             mag_array_1.SetValue(i, displacement_mag_1[i])

#         # Displacement magnitude from file 2
#         mag_array_2 = vtk.vtkFloatArray()
#         mag_array_2.SetName(f"{filename2}_displacement_magnitude")
#         mag_array_2.SetNumberOfValues(num_points)
#         for i in range(num_points):
#             mag_array_2.SetValue(i, displacement_mag_2[i])

#         # Absolute difference in displacement magnitudes
#         diff_array = vtk.vtkFloatArray()
#         diff_array.SetName("displacement_magnitude_difference")
#         diff_array.SetNumberOfValues(num_points)
#         for i in range(num_points):
#             diff_array.SetValue(i, displacement_mag_diff[i])

#         # Signed difference
#         signed_diff_vtk = vtk.vtkFloatArray()
#         signed_diff_vtk.SetName("displacement_magnitude_signed_difference")
#         signed_diff_vtk.SetNumberOfValues(num_points)
#         for i in range(num_points):
#             signed_diff_vtk.SetValue(i, signed_diff_array[i])

#         # Relative difference
#         relative_diff_vtk = vtk.vtkFloatArray()
#         relative_diff_vtk.SetName("displacement_magnitude_relative_difference")
#         relative_diff_vtk.SetNumberOfValues(num_points)
#         for i in range(num_points):
#             relative_diff_vtk.SetValue(i, relative_diff_array[i])

#         # Add arrays to output
#         output_pointData.AddArray(mag_array_1)
#         output_pointData.AddArray(mag_array_2)
#         output_pointData.AddArray(diff_array)
#         output_pointData.AddArray(signed_diff_vtk)
#         output_pointData.AddArray(relative_diff_vtk)

#         # Set the difference array as the active scalar for visualization
#         output_pointData.SetActiveScalars("displacement_magnitude_difference")
        
#         print(f"\nOutput arrays created:")
#         print(f"  - {filename1}_displacement_magnitude")
#         print(f"  - {filename2}_displacement_magnitude") 
#         print(f"  - displacement_magnitude_difference (ACTIVE SCALAR)")
#         print(f"  - displacement_magnitude_signed_difference")
#         print(f"  - displacement_magnitude_relative_difference")
#         print(f"  - Active scalar set to displacement_magnitude_difference for visualization")

#         # Write output file
#         writer = vtk.vtkUnstructuredGridWriter()
#         writer.SetInputData(output_data)
#         file_output = f"{filename1}_{filename2}_displacement_magnitude_diff.vtk"
#         output_full_filename = output_path + file_output
#         writer.SetFileName(output_full_filename)
#         writer.Update()
#         writer.Write()

#         # Calculate percentiles for outlier removal
#         p5 = np.percentile(displacement_mag_diff, 5)
#         p95 = np.percentile(displacement_mag_diff, 95)
        
#         # Filter out outliers (keep values between 5th and 95th percentiles)
#         filtered_diff = displacement_mag_diff[(displacement_mag_diff >= p5) & (displacement_mag_diff <= p95)]
        
#         # Enhanced statistics reporting
#         print(f"\nEnhanced displacement magnitude statistics:")
#         print(f"File 1 ({filename1}):")
#         print(f"  Min: {np.min(displacement_mag_1):.6f}")
#         print(f"  Max: {np.max(displacement_mag_1):.6f}")
#         print(f"  Mean: {np.mean(displacement_mag_1):.6f}")
#         print(f"  Std: {np.std(displacement_mag_1):.6f}")

#         print(f"File 2 ({filename2}):")
#         print(f"  Min: {np.min(displacement_mag_2):.6f}")
#         print(f"  Max: {np.max(displacement_mag_2):.6f}")
#         print(f"  Mean: {np.mean(displacement_mag_2):.6f}")
#         print(f"  Std: {np.std(displacement_mag_2):.6f}")

#         print(f"Signed Difference (File1 - File2):")
#         print(f"  Min: {np.min(signed_diff_array):.6f}")
#         print(f"  Max: {np.max(signed_diff_array):.6f}")
#         print(f"  Mean: {np.mean(signed_diff_array):.6f}")
#         print(f"  Std: {np.std(signed_diff_array):.6f}")

#         print(f"Absolute Difference:")
#         print(f"  Min: {np.min(displacement_mag_diff):.6f}")
#         print(f"  Max: {np.max(displacement_mag_diff):.6f}")
#         print(f"  Mean: {np.mean(displacement_mag_diff):.6f}")
#         print(f"  Std: {np.std(displacement_mag_diff):.6f}")

#         print(f"Relative Difference:")
#         print(f"  Min: {np.min(relative_diff_array):.6f}")
#         print(f"  Max: {np.max(relative_diff_array):.6f}")
#         print(f"  Mean: {np.mean(relative_diff_array):.6f}")
#         print(f"  Std: {np.std(relative_diff_array):.6f}")

#         # Check for systematic bias
#         positive_diffs = signed_diff_array[signed_diff_array > 0]
#         negative_diffs = signed_diff_array[signed_diff_array < 0]
#         zero_diffs = signed_diff_array[signed_diff_array == 0]

#         print(f"\nBias Analysis:")
#         print(f"  Points where File1 > File2: {len(positive_diffs)} ({100*len(positive_diffs)/num_points:.1f}%)")
#         print(f"  Points where File1 < File2: {len(negative_diffs)} ({100*len(negative_diffs)/num_points:.1f}%)")
#         print(f"  Points where File1 = File2: {len(zero_diffs)} ({100*len(zero_diffs)/num_points:.1f}%)")
#         if len(positive_diffs) > 0:
#             print(f"  Mean of positive differences: {np.mean(positive_diffs):.6f}")
#         if len(negative_diffs) > 0:
#             print(f"  Mean of negative differences: {np.mean(negative_diffs):.6f}")
        
#         print(f"Magnitude Difference (outliers removed - 5th to 95th percentiles):")
#         print(f"  5th percentile: {p5:.6f}")
#         print(f"  95th percentile: {p95:.6f}")
#         print(f"  Min (after outlier removal): {np.min(filtered_diff):.6f}")
#         print(f"  Max (after outlier removal): {np.max(filtered_diff):.6f}")
#         print(f"  Mean (after outlier removal): {np.mean(filtered_diff):.6f}")
#         print(f"  Number of points used: {len(filtered_diff)}/{len(displacement_mag_diff)} ({100*len(filtered_diff)/len(displacement_mag_diff):.1f}%)")

#         print(f"\nCOMPLETE: {output_full_filename}")
        
#     else:
#         print("Failed to open files, one or both of the following may not exist:")
#         print(fullFilename1)
#         print(fullFilename2)

def create_difference_plot(filename1, path_to_file_1, filename2, path_to_file_2, output_path, array_name=None):
    """
    Original function preserved for backward compatibility
    """
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

        assert numArray == numArray2, f"{numArray} and {numArray2} are not equal"
        assert data.GetNumberOfPoints() == data2.GetNumberOfPoints(), f"{data.GetNumberOfPoints()} and {data2.GetNumberOfPoints()} are not equal"

        # Create new unstructured grid
        output_data = vtk.vtkUnstructuredGrid()
        output_data.DeepCopy(data)
        output_pointData = output_data.GetPointData()
        output_pointData.Initialize()

        # Read in data and rename arrays
        for idx in range(0, numArray):
            arr_1 = data.GetPointData().GetArray(idx)
            arr_2 = data2.GetPointData().GetArray(idx)
            arr_full_name = arr_1.GetName()
            arr_full_name_2 = arr_2.GetName()
            
            if "material" not in arr_full_name:
                if array_name is None or arr_full_name in array_name:

                    assert arr_full_name == arr_full_name_2, f"{arr_full_name} and {arr_full_name_2} are not the same"

                    # Rename array from file 1
                    prop_name1 = filename1 + "_" + arr_full_name
                    
                    # Rename array from file 2
                    prop_name2 = filename2 + "_" + arr_full_name_2

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

                    # Create new arrays
                    ao = vtk.vtkFloatArray()
                    ao.SetName('Diff ' + arr_full_name)
                    ao.SetNumberOfValues(difference.size)
                    for x in range(difference.size):
                        ao.SetValue(x, difference[x])

                    ao1 = vtk.vtkFloatArray()
                    ao1.SetName(prop_name1)
                    ao1.SetNumberOfValues(arr_1_new.size)
                    for x in range(arr_1_new.size):
                        ao1.SetValue(x, arr_1_new[x])

                    ao2 = vtk.vtkFloatArray()
                    ao2.SetName(prop_name2)
                    ao2.SetNumberOfValues(arr_2_new.size)
                    for x in range(arr_2_new.size):
                        ao2.SetValue(x, arr_2_new[x])
                    
                    output_pointData.AddArray(ao)
                    output_pointData.AddArray(ao1)
                    output_pointData.AddArray(ao2)
                    print("Array " + ao.GetName() + " written")
            else:
                # Add material_id array
                output_pointData.AddArray(arr_1)
                print("Material array preserved")

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetInputData(output_data)
        file_output = filename1 + "_" + filename2 + ".vtk"
        output_full_filename = output_path + file_output
        writer.SetFileName(output_full_filename)
        writer.Update()
        writer.Write()

        print("COMPLETE: " + output_full_filename)
    else:
        print("Failed to open files, one or both of the following may not exist:")
        print(fullFilename1)
        print(fullFilename2)


if __name__ == "__main__":
    # Path configuration
    path = "./output_atrophy/post_processing/rampp_FA_9R/"
    out_path = path
    
    # Define all comparison pairs
    comparison_pairs = [
        ("rampp_atrophy_9R_last_timestep", "rampp_atrophy_FA_last_timestep")
        # ("OAS1_0004_MR1_atrophy_17R_last_timestep", "OAS1_0004_MR1_atrophy_2R_last_timestep"),
        # ("OAS1_0004_MR1_atrophy_17R_last_timestep", "OAS1_0004_MR1_atrophy_4R_last_timestep"),
        # ("OAS1_0004_MR1_atrophy_17R_last_timestep", "OAS1_0004_MR1_atrophy_9R_last_timestep")
    ]
    
    print("="*80)
    print("DISPLACEMENT MAGNITUDE DIFFERENCE ANALYSIS")
    print("="*80)
    
    # Process each comparison pair
    for i, (file1, file2) in enumerate(comparison_pairs, 1):
        print(f"\n{'='*20} COMPARISON {i}/4 {'='*20}")
        print(f"Comparing: {file1} vs {file2}")
        print("-" * 60)
        
        # Use the new function for displacement magnitude differences
        create_displacement_magnitude_difference_plot(file1, path, file2, path, out_path, 
                                                    displacement_array_name="displacement")
        
        if i < len(comparison_pairs):
            print("\n" + "="*60)
    
    print(f"\n{'='*80}")
    print("ALL COMPARISONS COMPLETED")
    print(f"{'='*80}")
    print(f"Output directory: {out_path}")
    print("Generated files:")
    for file1, file2 in comparison_pairs:
        file1_short = file1.split("_")[-2] + "_" + file1.split("_")[-1]
        file2_short = file2.split("_")[-2] + "_" + file2.split("_")[-1]
        print(f"  - {file1_short}_{file2_short}_displacement_magnitude_diff.vtk")