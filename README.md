## Description 17-region-model2 branch ##
This is an exension of the brain creation code, branch 17-region-model which was created by Nicole Tueni. The code within the branch ‘master’ does not work properly, so basically the 17-region model.
In addition to the atrophy code of the 17-region-model branch, this branch allows the processing of FA (fractional anisotropy) values within the meshing process. This mesh can then be used as an Input for the master atrophy branch of the EFI_atrophy_and_tumor branch.
To run the simulation, an aligned T1 image, a FA image and a aseg.mgz which can be retrieved by running the “recon -all” command in freesurfer is needed. The raw data coming from the MRI should be the “T1_iso.nii” (not the “t1_mprage_sag_p2” because the voxels might be anisotropic) and the “ep2d_diff_mddw30_p2_s2_2mm_FA.nii”. In order to align the MRI images, past them into one folder and paste the bash script “brain_transformation1” and the python file “brain_mask1” from the folder “Brain_Creation_Code-17-region-model\freesurfer_preprocessing\tools” to the same directory. Run the bash script (modify the working directory in line 9 and adjust the input file names in line 16-18) and “brain_mask1.py” afterwards to get the file
“aseg1.mgz”
In order to run the Brain_creation Code, adjust the project directory in the 
“Brain_Creation_Code-17-region-model\personal_scripts\create_model.py” line 5 and adjust the name of the output folder in line 31.
Finally, run “Brain_Creation_Code-17-region-model\personal_scripts\create_model.py”.


This project was started by Emma Griffiths ca. 2023 and continued by Nicole Tueni (see branches).

The code to create a model is found in the file 'create_model' in the 'personal_scripts' folder. Uncomment sections '####### ATROPHY CODE #######' and '####### TUMOR LESION CODE #######' as necessary. This code will both create a mesh from the MRI images and create the prm files according to the heterogeneity level specified. In creating the prm files the volume average of the mesh is used to create homogenized material properties.

To run a new brain creation code you need the following files: 1. The aparc segmented file from FreeSurfer. Please change the file name to 'aparc_DKTatlas+aseg.mgz' 2. The brain stem segmented file from Freesurfer. Please change the file name to 'brainstemSsLabels.mgz'

You will also need to change the location of the input files 'path_in' in lines 17 and 57 in "create_model". All files will be created in the a folder with their name found in the '/IOput/out/atrophy_files' or '/IOput/out/tumor_files' folder.

## Modifications ##
- Run python3 -i -m personal_scripts.create_model in Brain_Creation_Code directory
- Fixed problem with homogenization function
- Possible to get 1R - 2R - 4R - 9R and 17R models with updated material parameters
- Filled ventricles with fluid - CSF properties added to ventricles
- Smoothed out CSF
- Added BC option: Local boundary - allows to apply BC on one or more material labels
- Un-coarsen mesh: Preprocessor.py - super().coarsen()

