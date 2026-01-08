# Module Import

import sys
import os
import brain_creation
from config.Config import ConfigFile
import writers.HeterogeneityConverter as heterogConverter
from readers.Readers import Reader
from personal_scripts.create_prms import CreateAtrophyPRM, CreateTumorPRM
from personal_scripts.ucd_processing import process_UCD_extension 
from personal_scripts.ucd_processing import process_UCD_column_remove
from personal_scripts.ucd_processing import convert_to_vtk

# Model type options: basic_fullcsf, basic_partilacsf, basic_nocsf, atrophy, lesion

path_to_oasis = "./IOput/in/brain_FA/"
file_name_in = "aparc_DKTatlas+aseg.mgz"
path_to_out = "./IOput/out/atrophy_files/brain_FA"

filenames = ["brain_FA"] # change subject name here
for name in filenames:
    path_in = path_to_oasis
    path_out = "/".join([path_to_out, name])
    file_name_out = name
    
    config = ConfigFile(path_in, file_name_in, path_out, file_name_out,
                        configFilePath="./IOput/model_config.ini", model_type='atrophy')
    mesh = brain_creation.run(config)
    mesh.write(path_out, file_name_out, ['vtk','ucd'])
    # UNCOMMENT THE FOUR LINES BELOW AND COMMENT OUT LINES 21-24 inkl. IF YOU HAVE AN INPUT FIL ETO READ AND DON'T WANT TO CREATE A NEW MESH
    # reader = Reader('vtk')
    # reader.openReader(file_name_out + "_VTK", path_out)
    # mesh = reader.getMesh()
    # reader.closeReader()

    conditioning = 'preconditioned'
    poissons = '0,45'
    for heterogeneity_model in [heterogConverter.Heterogeneity.ONER, heterogConverter.Heterogeneity.TWOR, heterogConverter.Heterogeneity.FOURR,
                                heterogConverter.Heterogeneity.NINER, heterogConverter.Heterogeneity.NINETEENR]:
    # for heterogeneity_model in [heterogConverter.Heterogeneity.NINETEENR]:

        atrophy_creator = CreateAtrophyPRM("./personal_scripts/atrophy_template_folder/atrophy_template_V2.prm")
        atrophy_creator.create_materials(mesh, conditioning, poissons, heterogeneity_model)
        atrophy_creator.write_materials()
        atrophy_creator.complete_prm(path_out, file_name_out, "{}_atrophy_{}R".format(file_name_out, heterogeneity_model.value))
        output_prm = "/".join([path_out, "{}_atrophy_{}R".format(file_name_out, heterogeneity_model.value)])
        atrophy_creator.write_prm(output_prm)
        atrophy_creator.close_prm()

    # change add and remove columns for the final output files
    input_file_extension = os.path.join(path_out, "brain_FA_UCD.inp")
    output_file_extension = os.path.join(path_out, "brain_FA_UCD1.inp")
    process_UCD_extension(input_file_extension, output_file_extension)


    input_file_column = os.path.join(path_out, "brain_FA_UCD.inp")
    output_file_column = os.path.join(path_out, "brain_FA_UCD_orig.inp")
    process_UCD_column_remove(input_file_column, output_file_column)

    input_file_column = os.path.join(path_out, "brain_FA_UCD1.inp")
    output_file_column = os.path.join(path_out, "brain_FA_UCD_FA.inp")
    process_UCD_column_remove(input_file_column, output_file_column)

    input_file = os.path.join(path_out, "brain_FA_UCD_FA.inp")
    output_file = os.path.join(path_out, "brain_FA_VTK_FA.vtk")
    convert_to_vtk(input_file, output_file)

    print("COMPLETE")
    print("Files written to",path_out)

