import numpy as np
import pandas as pd
import vtk
import os
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # headless
import matplotlib.pyplot as plt
from os.path import exists

# --- CONFIGURATION ---
base_path = "./output_atrophy/post_processing/rampp_final/"
fileNames = [
    "rampp_atrophy_9R_last_timestep",
    "rampp_atrophy_FA_last_timestep"
]
hist_output_dir = os.path.join(base_path, "histograms_all_measures")
os.makedirs(hist_output_dir, exist_ok=True)
measure = [
    'max_principal_stretch', 'med_principal_stretch', 'min_principal_stretch', 'max_shear',  'hyrdostatic_stress', "von_mises"]
metric_labels = {
    'max_principal_stretch': 'First principal stretch [-]',
    'med_principal_stretch': 'Second principal stretch [-]',
    'min_principal_stretch': 'Third principal stretch [-]',
    'max_shear': 'Maximum shear [-]',
    'hyrdostatic_stress': 'Hydrostatic stresses [MPa]',
    'von_mises': 'von Mises stresses [MPa]',
    'displacement': 'Displacement'
}

regions = {
    'Corpus Callosum': [251],
    'Corona Radiata': [2],
    'Hippocampus': [17],
    'Cerebellum': [7],
    'Midbrain':  [173, 10],
    'Brain Stem': [174, 175],
    'Basal Ganglia': [11, 13, 12],
    'Amygdala': [18],
    'Cortex':  [1030, 1035, 1028, 1011, 1024],
    'Ventricles': [4]  # excluded from exports/plots
}

region_id_to_label = {v[0]: k for k, v in regions.items()}

# --- X-AXIS CROPPING CONFIGURATION ---
CROP_METHOD = 'percentile'  # Options: 'percentile', 'density_threshold', 'data_percentage'
PERCENTILE_RANGE = (1, 99)  # Crop to this percentile range
# DENSITY_THRESHOLD = 0.01  # For density method: minimum density as fraction of max
# DATA_PERCENTAGE = 0.98  # For data percentage method: include this fraction of data

def find_crop_limits(data1, data2, method='percentile'):
    """Find x-axis limits to crop out low-frequency noise"""
    combined_data = np.concatenate([data1, data2])
    
    if method == 'percentile':
        lower_limit = np.percentile(combined_data, PERCENTILE_RANGE[0])
        upper_limit = np.percentile(combined_data, PERCENTILE_RANGE[1])
        
    elif method == 'density_threshold':
        # Create histogram and find where density drops below threshold
        hist, bin_edges = np.histogram(combined_data, bins=200, density=True)
        max_density = np.max(hist)
        threshold = max_density * DENSITY_THRESHOLD
        
        # Find first and last bins above threshold
        above_threshold = hist >= threshold
        if np.any(above_threshold):
            first_idx = np.where(above_threshold)[0][0]
            last_idx = np.where(above_threshold)[0][-1]
            lower_limit = bin_edges[first_idx]
            upper_limit = bin_edges[last_idx + 1]
        else:
            lower_limit = np.min(combined_data)
            upper_limit = np.max(combined_data)
            
    elif method == 'data_percentage':
        # Include central percentage of data
        exclude_percentage = (1 - DATA_PERCENTAGE) / 2 * 100
        lower_limit = np.percentile(combined_data, exclude_percentage)
        upper_limit = np.percentile(combined_data, 100 - exclude_percentage)
    
    else:
        lower_limit = np.min(combined_data)
        upper_limit = np.max(combined_data)
    
    return lower_limit, upper_limit

def load_vtk_data(fullFilename):
    if not exists(fullFilename):
        print(f"File {fullFilename} does not exist. Skipping.")
        return None
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(fullFilename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    num_points = data.GetNumberOfPoints()
    num_arrays = data.GetPointData().GetNumberOfArrays()
    material_ids = None
    filenameData = {}
    for idx in range(num_arrays):
        arr_1 = data.GetPointData().GetArray(idx)
        arr_full_name = arr_1.GetName()
        if "material" not in arr_full_name:
            py_arr_1 = [arr_1.GetValue(i) for i in range(num_points)]
            filenameData[arr_full_name] = py_arr_1
        else:
            material_ids = [arr_1.GetValue(i) for i in range(num_points)]
    if material_ids is None:
        print(f"No material IDs found in {fullFilename}. Skipping.")
        return None
    region_data = {}
    for key, array in filenameData.items():
        region_values = {}
        for i, mat_id in enumerate(material_ids):
            region_values.setdefault(mat_id, []).append(array[i])
        region_data[key] = region_values
    return region_data

# --- DATA LOADING AND ORGANIZATION ---
all_data = []
for f in fileNames:
    full_path = os.path.join(base_path, f + ".vtk")
    region_data = load_vtk_data(full_path)
    if region_data is None:
        continue
    model_name = f.split('_')[-3]
    for metric in measure:
        if metric not in region_data:
            continue
        for region_id, values in region_data[metric].items():
            region_label = region_id_to_label.get(region_id)
            if region_label is None:
                continue
            for v in values:
                all_data.append({
                    'Region': region_label,
                    'Value': v,
                    'Simulation': model_name,
                    'Metric': metric
                })

# --- HISTOGRAM PLOTTING WITH X-AXIS CROPPING ---
region_order = list(regions.keys())

for selected_metric in measure:
    df = pd.DataFrame([d for d in all_data if d['Metric'] == selected_metric])
    if df.empty:
        print(f"No data for metric: {selected_metric}")
        continue

    for region in region_order:
        if region not in df['Region'].values:
            print(f"Warning: Region {region} not found in data")
            continue

        v1 = df.loc[(df.Region == region) & (df.Simulation == '9R'), 'Value'].values
        v2 = df.loc[(df.Region == region) & (df.Simulation == 'FA'), 'Value'].values

        if len(v1) == 0 or len(v2) == 0:
            print(f"Warning: No data found for region {region}")
            continue

        # Remove NaN values
        v1 = v1[~np.isnan(v1)]
        v2 = v2[~np.isnan(v2)]

        if len(v1) == 0 or len(v2) == 0:
            print(f"Warning: No valid data after NaN removal for region {region}")
            continue

        # Find crop limits to exclude low-frequency noise
        x_min, x_max = find_crop_limits(v1, v2, method=CROP_METHOD)
        
        # Calculate means (using all data, not just cropped)
        mean_v1 = np.nanmean(v1)
        mean_v2 = np.nanmean(v2)

        print(f"\nProcessing region {region}: 9R={len(v1)} points, FA={len(v2)} points")
        print(f"  X-axis cropped from [{np.min(np.concatenate([v1, v2])):.3f}, {np.max(np.concatenate([v1, v2])):.3f}] to [{x_min:.3f}, {x_max:.3f}]")

        # Create histogram (overlaid) with cropped x-axis
        plt.figure(figsize=(8, 5))
        
        # Create bins within the cropped range
        bins = np.linspace(x_min, x_max, 500)
        
        plt.hist(v1, bins=bins, alpha=0.7, label=f'Region-based (mean={mean_v1:.3f})', color=(230/255, 97/255, 0/255))
        plt.hist(v2, bins=bins, alpha=0.7, label=f'Voxel-based (mean={mean_v2:.3f})', color=(93/255, 58/255, 155/255))
        
        # Only show mean lines if they're within the cropped range
        if x_min <= mean_v1 <= x_max:
            plt.axvline(mean_v1, color=(33/255, 113/255, 181/255), linestyle='dashed', linewidth=2)
        if x_min <= mean_v2 <= x_max:
            plt.axvline(mean_v2, color=(255/255, 255/255, 0), linestyle='dashed', linewidth=2)
        
        # Set x-axis limits to crop out noise
        plt.xlim(x_min, x_max)
        
        plt.legend(loc='upper left', fontsize=14)
        plt.xlabel(metric_labels.get(selected_metric, selected_metric), fontsize=20, fontweight='bold')
        plt.ylabel("Frequency", fontsize=20, fontweight='bold')

        # Annotate means in a box
        mean_text = (f"Mean 9R = {mean_v1:.3f}\n"
                     f"Mean FA = {mean_v2:.3f}")
        plt.gca().text(0.98, 0.98, mean_text,
                       transform=plt.gca().transAxes,
                       verticalalignment='top',
                       horizontalalignment='right',
                       bbox=dict(boxstyle='round', facecolor='none', alpha=0.8),
                       fontsize=11)

        filename = os.path.join(hist_output_dir, f"hist_{selected_metric}_{region}_cropped.png")
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved cropped histogram: {filename}")

        # Difference histogram (v2 - v1, pairwise if possible) with cropping
        min_len = min(len(v1), len(v2))
        if min_len > 0:
            diff = v2[:min_len] - v1[:min_len]
            diff = diff[~np.isnan(diff)]  # Remove NaN values
            
            if len(diff) > 0:
                mean_diff = np.nanmean(diff)
                
                # Find crop limits for difference data
                diff_min = np.percentile(diff, PERCENTILE_RANGE[0])
                diff_max = np.percentile(diff, PERCENTILE_RANGE[1])
                
                plt.figure(figsize=(8, 5))
                bins_diff = np.linspace(diff_min, diff_max, 500)
                plt.hist(diff, bins=bins_diff, alpha=0.7, color='gray', 
                        label=f'Difference (mean={mean_diff:.3f})')
                
                # Only show mean line if it's within the cropped range
                if diff_min <= mean_diff <= diff_max:
                    plt.axvline(mean_diff, color='black', linestyle='dashed', linewidth=2)
                
                # Set x-axis limits to crop out noise
                plt.xlim(diff_min, diff_max)
                
                plt.legend(loc='upper left')
                plt.xlabel("Difference", fontsize=20, fontweight='bold')
                plt.ylabel("Frequency", fontsize=20, fontweight='bold')

                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                
                # Annotate mean difference in a box
                mean_diff_text = f"Mean diff = {mean_diff:.3f}"
                plt.gca().text(0.98, 0.98, mean_diff_text,
                               transform=plt.gca().transAxes,
                               verticalalignment='top',
                               horizontalalignment='right',
                               bbox=dict(boxstyle='round', facecolor='none', alpha=0.8),
                               fontsize=11)
                
                diff_filename = os.path.join(hist_output_dir, f"hist_diff_{selected_metric}_{region}_cropped.png")
                plt.tight_layout()
                plt.savefig(diff_filename, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"Saved cropped difference histogram: {diff_filename}")

print(f"\nX-axis cropping configuration used:")
print(f"- Crop method: {CROP_METHOD}")
if CROP_METHOD == 'percentile':
    print(f"- Percentile range: {PERCENTILE_RANGE}")
elif CROP_METHOD == 'density_threshold':
    print(f"- Density threshold: {DENSITY_THRESHOLD}")
elif CROP_METHOD == 'data_percentage':
    print(f"- Data percentage included: {DATA_PERCENTAGE}")