import os
import nibabel as nib
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))

fa_img_path = os.path.join(script_dir, 'output_resampled.mgz')
fa_img = nib.load(fa_img_path)
mask_path = os.path.join(script_dir, 'brain_mask.mgz')
mask_img = nib.load(mask_path)


# 2. Convert the image data to Numpy arrays
fa_data = fa_img.get_fdata()  # FA-array
mask_data = mask_img.get_fdata()  # brain amsk array

# 3. Use brain mask to put all points outside the brain to 0
masked_fa_data = fa_data * mask_data  # only the values inside the mask stay


#4. Create new NIFTI image with masked FA-values
masked_fa_img = nib.Nifti1Image(masked_fa_data, fa_img.affine, fa_img.header)

# 5. Save the masked image
output_path = os.path.join(script_dir, 'aseg1.mgz')
nib.save(masked_fa_img, output_path)

print("Masking finished. Output file: 'aseg1.mgz'.")
