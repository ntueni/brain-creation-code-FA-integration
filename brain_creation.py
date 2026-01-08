# -*- coding: utf-8 -*-
"""
Brain creator 

Created on Thu Jun 15 09:34:00 2023

@author: grife

This script runs the full execution required to create a 3D brain model from a
freesurfer aseg file.

Input file given by path + fileIn
Preferences for model creation are defined in ConfigFile class
Output is a folder Models created in current working directory

main - the main function of the script
"""
import save_data
import numpy as np
from BrainHexMesh import BrainHexMesh
from mesh.PostProcessorFactory import PostProcessorFactory
from point_cloud.PointCloud import PointCloud
from config.Config import ConfigFile
from voxel_data import Preprocessor
from writers.aseg_manipulate import create_aseg
from writers.Writers import Writer
from ArrayProcessor import ArrayProcessor

def run(config):

    brainCreator = BrainHexMesh(config)

    # Writes configuration preferences to output location
    config.open_config_file()
    config.write_preamble()

    # Gets aseg data as 3D of segmentation labels in voxels
    data = brainCreator.import_data()

    # Homogenizes data label according to Materials label as defined in Config class
    data = brainCreator.homogenize_data(data)

    # Pre-processes data to ensure valid mesh base on config setting:
    preprocessor = Preprocessor.PreProcessorFactory.get_preprocessor(config, data)
    assert isinstance(preprocessor, Preprocessor.IPreprocessor)
    data_new = preprocessor.preprocess_data()
    create_aseg(config.get("file_in_path"), config.get("file_in"), config.get("file_out_path"), config.get("fileout"),
                data_new)

   # Creates point cloud from voxel data

    """"
    newData_4d einfügen und mit newData vereinen.
    """
    newData_4d = save_data.get_data()
    if newData_4d is not None:
        print(f"newData_4d successfully loaded with shape: {newData_4d.shape}")
    else:
        print("No data available in newData_4d. Ensure it was saved before running this.")

    # Neues 4D-Array initialisieren
    combined_array = np.zeros(data_new.shape[:3] + (2,))
    
    # Materialwerte in den ersten Kanal kopieren
    combined_array[..., 0] = data_new[..., -1]  
    
    # FA-Werte nur für übereinstimmende Koordinaten in den zweiten Kanal übernehmen
    for x in range(data_new.shape[0]):
        for y in range(data_new.shape[1]):
            for z in range(data_new.shape[2]):
                if data_new[x, y, z] != 0:  # Materialwert ist nicht null
                    combined_array[x, y, z, 0] = data_new[x, y, z]
                    if newData_4d[x, y, z, 1] != 0:  # FA-Wert ist nicht null
                        combined_array[x, y, z, 1] = newData_4d[x, y, z, 1]
                    else:  # FA-Wert ist null
                        combined_array[x, y, z, 1] = 0
                else:  # Materialwert ist null
                    combined_array[x, y, z, 1] = 0  # FA-Wert ebenfalls auf null setzen
    
    processor = ArrayProcessor()
    processor.printArray(combined_array, "combined_array")
    data_new = combined_array
    
    print(f"shape data_new: {data_new.shape}" )

    pointCloud = PointCloud()
    pointCloud.create_point_cloud_from_voxel(data_new)

    # Creates mesh from point cloud
    mesh = brainCreator.make_mesh(pointCloud.pcd)
    brainCreator.clean_mesh(mesh)

    # Moves mesh to the center of the corpus callosum
    mesh.center_mesh_by_region(251)

    # Wrapping of post-processing operations (operation selection defined in config file)
    post_processor = PostProcessorFactory.get_post_processor(mesh, config)
    post_processor.post_process()

    # Close config file write out
    config.close_config_file()

    return mesh


if __name__ == "__main__":
    # Model type options: basic_fullcsf, basic_partilacsf, basic_nocsf, atrophy, lesion
    config = ConfigFile("./IOput/in", "aseg.mgz","./IOput/out",
                        "brain_smoothed", model_type='basic_nocsf')
    mesh = run(config)
    mesh.write(config.get('file_out_path'), config.get('fileout'))



