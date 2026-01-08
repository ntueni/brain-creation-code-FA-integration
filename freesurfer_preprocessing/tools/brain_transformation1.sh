#!/bin/bash
# NIfTI files are converted into the mgz format, registered, resampled and then provided with a brain mask.
# Prerequisites: FreeSurfer and Python (with nibabel and numpy) must be installed.
# The input files (T1.nii, ep2d_diff_mddw_20_p2_FA.nii and aseg.mgz) must be in the current working directory.
# after running the script you have to run brain_mask.py file

set -e
# adjust working directory here
WORK_DIR="/mnt/c/Users/pumab/Documents/6_Brain_segmentation/rampp/brain_segmentation2"

cd "$WORK_DIR" || { echo "Error: Could not change to working directory."; exit 1; }

echo "Working directory: $(pwd)"

# files needed
T1_INPUT="T1_iso.nii"
FA_INPUT="ep2d_diff_mddw_20_p2_FA.nii"
ASEG_INPUT="aseg.mgz"

echo "Eingangsdateien:"
echo "T1: $T1_INPUT"
echo "FA: $FA_INPUT"
echo "aseg: $ASEG_INPUT"

echo "--------------------------------------------------"
echo "Step 1.1: Convert and register the T1 file"
echo "--------------------------------------------------"
# Convert NIfTI-Dateien to the mgz format
mri_convert "$T1_INPUT" T1.mgz
mri_convert "$FA_INPUT" FA_converted.mgz

echo "--------------------------------------------------"
echo "Step 1.2: Align the T1 and the FA-picture"
echo "--------------------------------------------------"

#register FA to the T1 picture (create FA_aligned.mgz and FA_aligned.lta)
mri_robust_register --mov FA_converted.mgz --dst T1.mgz --lta FA_aligned.lta --mapmov FA_aligned.mgz --satit --maxit 10

echo "--------------------------------------------------"
echo "Step 1.3: Resample of the registered FA file"
echo "--------------------------------------------------"
# Resample of the registered FA file to the same resolution as aseg.mgz (256x256x256)
mri_convert --resample_type cubic --like "$ASEG_INPUT" FA_aligned.mgz output_resampled.mgz

echo "--------------------------------------------------"
echo "Step 2: Create a brain mask"
echo "--------------------------------------------------"
# Create a brain mask from aseg.mgz (all voxels with a value >= 1)
mri_binarize --i "$ASEG_INPUT" --o brain_mask.mgz --min 1

echo "--------------------------------------------------"
echo "Step 2.1: Run the Python file brain_mask.py in order to get the final FA-File aseg1.mgz"
echo "--------------------------------------------------"
