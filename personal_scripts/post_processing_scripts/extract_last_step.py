import os
import xml.etree.ElementTree as ET
import vtk


### This code allows to extract the last step from the .pvd file
### containing the atrophy/tumor results.
### Returns vtk files filename.vtk


# Set the paths
base_directory = "./output_atrophy/rampp_final"
output_path =  "./output_atrophy/post_processing/rampp_final/"

for folder_name in os.listdir(base_directory):
    folder_path = os.path.join(base_directory, folder_name)
    
    # Change folder name in line below 
    if os.path.isdir(folder_path) and "rampp" in folder_name:
        # Construct the path to the .pvd file inside the current folder
        pvd_file_name = f"{folder_name}_solution-3d.pvd"
        print(pvd_file_name)
        pvd_file_path = os.path.join(folder_path, pvd_file_name)
        print(pvd_file_path)

        if os.path.isfile(pvd_file_path):
            # Parse the .pvd file to extract the last timestep
            print(f"Processing: {pvd_file_path}")
            tree = ET.parse(pvd_file_path)
            root = tree.getroot()
            
            # Get the list of all <DataSet> elements under <Collection>
            datasets = root.findall(".//DataSet")
            if not datasets:
                print(f"No datasets found in {pvd_file_path}. Skipping...")
                continue
            
            # Extract the last dataset (assumed to be the last timestep)
            last_dataset = datasets[-1]
            vtu_filename = last_dataset.get("file")
            vtu_file_path = os.path.join(folder_path, vtu_filename)
            
            if os.path.isfile(vtu_file_path):
                print(f"Extracting last timestep from: {vtu_file_path}")
                
                # Read the .vtu file using a generic VTK reader
                reader = vtk.vtkXMLGenericDataObjectReader()
                reader.SetFileName(vtu_file_path)
                reader.Update()
                
                # Get the output data
                output_data = reader.GetOutput()
                
                # Debug: Check if the output data has any content
                num_points = output_data.GetNumberOfPoints() if hasattr(output_data, "GetNumberOfPoints") else 0
                num_cells = output_data.GetNumberOfCells() if hasattr(output_data, "GetNumberOfCells") else 0
                print(f"Number of points in the last timestep: {num_points}")
                print(f"Number of cells in the last timestep: {num_cells}")
                
                if num_points == 0 or num_cells == 0:
                    print(f"Warning: Empty or unsupported data in {vtu_file_path}. Skipping...")
                    continue
                
                # Save the output to a .vtk file if data is present
                vtk_output_file = os.path.join(output_path, f"{folder_name}_last_timestep.vtk")
                writer = vtk.vtkDataSetWriter()
                writer.SetFileName(vtk_output_file)
                writer.SetInputData(output_data)
                writer.Write()
                
                print(f"Saved last timestep to {vtk_output_file}")
            else:
                print(f"VTU file not found: {vtu_file_path}")
