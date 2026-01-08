# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Generate volume-change plots for BOTH mappings:
# 1) Original REGION_NAME_MAPPING (per-region)
# 2) Grouped mapping you requested (Am, BG, BS, C, Hi, CC, CR, M, CB)

# Usage:
#     python volume_change_plots_both_mappings.py

# Notes:
# - Expects VTK files under INPUT_DIRECTORY_PATH with filenames: {BASE}_{step}.vtk
#   for BASE in BASE_FILENAMES and step in 0..16 (as in your original script).
# - Outputs PNG plots into OUTPUT_DIRECTORY_PATH.
# """

# import os
# import sys
# import vtk
# import numpy as np
# from vtk.util.numpy_support import vtk_to_numpy
# import matplotlib
# matplotlib.use('Agg')  # headless
# import matplotlib.pyplot as plt

# # -------------------- Config --------------------
# INPUT_DIRECTORY_PATH = './output_atrophy/post_processing/rampp_FA_9R/brainrampp/'
# OUTPUT_DIRECTORY_PATH = './output_atrophy/post_processing/rampp/plots/volume_evolution_per_region/'

# BASE_FILENAMES = ['rampp9R', 'ramppFA']  # models
# TIMESTEPS = list(range(17))  # 0..16

# # Original mapping setup
# BRAIN_REGIONS = [1011, 7, 2, 174, 173, 175, 1030, 17, 1024, 1028, 11, 251, 13, 12, 1035, 10, 18, 4]
# REGION_NAME_MAPPING = {
#     1024: "MC", 1030: "TL", 1028: "FC", 2: "CR", 1035: "CI", 12: "Pu", 13: "Pa", 10: "Th",
#     251: "CC", 11: "NC", 18: "Am", 17: "Hi", 173: "M", 174: "P", 1011: "VC", 7: "CB", 175: "Me", 4: "V"
# }


# GROUPED_REGIONS = {
#     'Corpus Callosum': [251],               # CC = CC
#     'Corona Radiata': [2],                 # CR = CR
#     'Hippocampus': [17],                # Hi = Hi
#     'Cerebellum': [7],                  # CB = CB
#     'Midbrain':  [173, 10],           # M = M, Th
#     'Brain Stem': [174, 175],          # BS = P, Me
#     'Basal Ganglia': [11, 13, 12],        # BG = NC, Pa, Pu   
#     'Amygdala': [18],                # Am = Am
#     'Cortex':  [1030, 1035, 1028, 1011, 1024],  # C = TL, CI, FC, VC, MC
#     'Ventricles': [4]                  # V = V
# }

# # Colors (Oxford blue & OI yellow)
# COLOR_9R = np.array([0, 114, 178]) / 255.0
# COLOR_FA = np.array([240, 228, 66]) / 255.0

# os.makedirs(OUTPUT_DIRECTORY_PATH, exist_ok=True)

# # -------------------- Helpers --------------------

# def compute_hex_volume(coords: np.ndarray) -> float:
#     """Approximate hexahedral cell volume from 8 corner points (as in your script)."""
#     if len(coords) == 8:
#         v0 = coords[0]
#         return abs(np.dot(np.cross(coords[1] - v0, coords[2] - v0), coords[4] - v0)) / 6.0
#     return 0.0

# def read_material_array(ugrid) -> np.ndarray:
#     """Return the material id array (point data) as numpy array; raise if missing."""
#     material_ids_array = None
#     pd = ugrid.GetPointData()
#     for idx in range(pd.GetNumberOfArrays()):
#         arr = pd.GetArray(idx)
#         if arr and arr.GetName() and "material" in arr.GetName().lower():
#             material_ids_array = arr
#             break
#     if material_ids_array is None:
#         raise RuntimeError("Material IDs array not found in VTK file.")
#     return vtk_to_numpy(material_ids_array)

# def load_ugrid(vtk_path: str):
#     """Read a legacy VTK unstructured grid file."""
#     reader = vtk.vtkUnstructuredGridReader()
#     reader.SetFileName(vtk_path)
#     reader.Update()
#     return reader.GetOutput()

# def accumulate_volumes(ugrid, material_ids_np, region_ids):
#     """Compute total volume per material id in region_ids and total/24 bookkeeping."""
#     volumes = {rid: 0.0 for rid in region_ids}
#     volumes.update({'entire_volume': 0.0, 24: 0.0})
#     n_cells = ugrid.GetNumberOfCells()
#     for cid in range(n_cells):
#         cell = ugrid.GetCell(cid)
#         pid0 = cell.GetPointIds().GetId(0)
#         mat_id = int(material_ids_np[pid0])
#         pts = cell.GetPoints()
#         coords = np.array([pts.GetPoint(i) for i in range(pts.GetNumberOfPoints())], dtype=float)
#         v = compute_hex_volume(coords)
#         volumes['entire_volume'] += v
#         if mat_id in volumes:
#             volumes[mat_id] += v
#         # 24 bookkeeping
#         if 24 in volumes and mat_id == 24:
#             volumes[24] += v
#     return volumes

# def summarize_last_step_changes(all_data, label_map, outfile_prefix, title_suffix=""):
#     """Make grouped bar plots (unsorted and sorted) for last-step percentage changes."""
#     # Extract last-step changes & labels
#     labels = []
#     changes_9R = []
#     changes_FA = []
#     for key in label_map:  # keep the natural order of label_map
#         labels.append(key)
#         changes_9R.append(all_data['rampp9R'][key][-1])
#         changes_FA.append(all_data['ramppFA'][key][-1])

#     # Helper to label bars
#     def add_value_labels(ax, bars, values):
#         for bar, value in zip(bars, values):
#             h = bar.get_height()
#             ax.text(bar.get_x() + bar.get_width()/2., 
#                     h + (0.5 if h >= 0 else -1.0),
#                     f'{value:.1f}%',
#                     ha='center', va='bottom' if h >= 0 else 'top', fontsize=9)

#     # --- Unsorted ---
#     x = np.arange(len(labels))
#     bw = 0.35
#     fig, ax = plt.subplots(figsize=(16, 10))
#     b1 = ax.bar(x - bw/2, changes_9R, bw,
#             color=COLOR_9R, edgecolor='black',
#             label='9R region based μ', alpha=0.85, linewidth=0.8)

#     b2 = ax.bar(x + bw/2, changes_FA, bw,
#             color=COLOR_FA, edgecolor='black',
#             label='FA voxel based μ', alpha=0.85, linewidth=0.8)
#     add_value_labels(ax, b1, changes_9R); add_value_labels(ax, b2, changes_FA)
#     # ax.set_xlabel("Brain Regions" if title_suffix == "Original" else "Brain Regions", fontsize=20, fontweight='bold')
#     # ax.set_ylabel("Percentage Change in Volume (%)", fontsize=20, fontweight='bold')
#     ax.set_ylabel("Percentage Change in Volume (%)", fontsize=20)

#     ax.set_xticks(x, labels, rotation=45, ha='right', fontsize=14)
#     plt.yticks(fontsize=14)
#     # ax.grid(True, which='both', axis='y', linestyle='--', alpha=0.6)
#     ax.legend(fontsize=18, loc='lower right')
#     fig.tight_layout()
#     out_unsorted = os.path.join(OUTPUT_DIRECTORY_PATH, f"{outfile_prefix}_Last_Step_Change.png")
#     fig.savefig(out_unsorted, dpi=300, bbox_inches='tight')
#     plt.close(fig)

#     # --- Sorted by average of 9R and FA ---
#     combo = list(zip(labels, changes_9R, changes_FA))
#     combo.sort(key=lambda t: (t[1] + t[2]) / 2.0)
#     slabels, s9, sFA = zip(*combo) if combo else ([], [], [])
#     x = np.arange(len(slabels))
#     fig, ax = plt.subplots(figsize=(16, 10))
#     b1 = ax.bar(x - bw/2, changes_9R, bw,
#                 color=COLOR_9R, edgecolor='black',
#                 label='9R region based μ', alpha=0.85, linewidth=0.8)

#     b2 = ax.bar(x + bw/2, changes_FA, bw,
#                 color=COLOR_FA, edgecolor='black',
#                 label='FA voxel based μ', alpha=0.85, linewidth=0.8)
#     add_value_labels(ax, b1, s9); add_value_labels(ax, b2, sFA)
#     ax.set_xlabel("Brain Regions" if title_suffix == "Original" else "Brain Regions", fontsize=20)
#     ax.set_ylabel("Percentage Change at Last Step (%)", fontsize=20)
#     ax.set_xticks(x, labels, rotation=45, ha='right', fontsize=14,)
#     plt.yticks(fontsize=14)
#     # ax.grid(True, which='both', axis='y', linestyle='--', alpha=0.2)
#     ax.legend(fontsize=18, loc='lower right')
#     fig.tight_layout()
#     out_sorted = os.path.join(OUTPUT_DIRECTORY_PATH, f"{outfile_prefix}_Last_Step_Change_Sorted.png")
#     fig.savefig(out_sorted, dpi=300, bbox_inches='tight')
#     plt.close(fig)

#     print(f"Saved: {out_unsorted}")
#     print(f"Saved: {out_sorted}")

# # -------------------- Main processing --------------------
# def main():
#     # Structures to hold % changes per timestep for ORIGINAL region labels (by display name)
#     per_region = {base: {REGION_NAME_MAPPING[r]: [] for r in BRAIN_REGIONS} for base in BASE_FILENAMES}
#     # Structures to hold % changes per timestep for GROUPED labels
#     per_group = {base: {g: [] for g in GROUPED_REGIONS.keys()} for base in BASE_FILENAMES}

#     for base in BASE_FILENAMES:
#         print(f"\nProcessing model: {base}")
#         # initial absolute volumes for each ORIGINAL region / GROUP for step 0
#         init_region_vols = None
#         init_group_vols = None

#         for step in TIMESTEPS:
#             vtk_path = os.path.join(INPUT_DIRECTORY_PATH, f"{base}_{step}.vtk")
#             if not os.path.exists(vtk_path):
#                 print(f"  Missing file: {vtk_path} (skipping this step)")
#                 continue

#             print(f"  Step {step}: {vtk_path}")
#             ugrid = load_ugrid(vtk_path)
#             mat_np = read_material_array(ugrid)

#             # --- ORIGINAL region volumes at this step ---
#             vols_region = accumulate_volumes(ugrid, mat_np, BRAIN_REGIONS)
#             # Brain volume excluding 24 (if needed elsewhere); computed but not used further
#             _brain_volume = vols_region['entire_volume'] - vols_region.get(24, 0.0)

#             # Track absolute per-region
#             abs_region = {REGION_NAME_MAPPING[r]: vols_region.get(r, 0.0) for r in BRAIN_REGIONS}

#             # --- GROUPED volumes at this step ---
#             # compute once by summing member material IDs
#             abs_group = {}
#             for gname, members in GROUPED_REGIONS.items():
#                 abs_group[gname] = float(sum(vols_region.get(mid, 0.0) for mid in members))

#             # Initialize baselines at first encountered step
#             if init_region_vols is None:
#                 init_region_vols = abs_region.copy()
#             if init_group_vols is None:
#                 init_group_vols = abs_group.copy()

#             # Percentage changes relative to step-0
#             for lbl in abs_region:
#                 base_v = init_region_vols[lbl]
#                 delta_pct = 0.0 if base_v == 0 else 100.0 * (abs_region[lbl] - base_v) / base_v
#                 per_region[base][lbl].append(delta_pct)

#             for g in abs_group:
#                 base_v = init_group_vols[g]
#                 delta_pct = 0.0 if base_v == 0 else 100.0 * (abs_group[g] - base_v) / base_v
#                 per_group[base][g].append(delta_pct)

#     # --- Plots for ORIGINAL mapping ---
#     print("\nGenerating plots for ORIGINAL mapping (last-step %) ...")
#     summarize_last_step_changes(per_region, 
#                                 label_map=[REGION_NAME_MAPPING[r] for r in BRAIN_REGIONS],
#                                 outfile_prefix="OriginalMapping",
#                                 title_suffix="Original")

#     # --- Plots for GROUPED mapping ---
#     print("\nGenerating plots for GROUPED mapping (last-step %) ...")
#     summarize_last_step_changes(per_group, 
#                                 label_map=list(GROUPED_REGIONS.keys()),
#                                 outfile_prefix="GroupedMapping",
#                                 title_suffix="Grouped")
    

#         # --- Export CSV tables for LaTeX/pgfplots (9 grouped regions, last-step %) ---
#     regions_order = list(GROUPED_REGIONS.keys())

#     def last_val(dct, base, key):
#         vals = dct[base][key]
#         return vals[-1] if vals else np.nan

#     # Long format: Region,Model,PercentChange
#     csv_long = os.path.join(OUTPUT_DIRECTORY_PATH, "GroupedMapping_LastStep_PercentChange_LONG.csv")
#     with open(csv_long, "w", encoding="utf-8") as f:
#         f.write("Region,Model,PercentChange\n")
#         for r in regions_order:
#             f.write(f"{r},9R,{last_val(per_group, 'rampp9R', r):.6f}\n")
#             f.write(f"{r},FA,{last_val(per_group, 'ramppFA', r):.6f}\n")

#     # Wide format: Region,9R,FA
#     csv_wide = os.path.join(OUTPUT_DIRECTORY_PATH, "GroupedMapping_LastStep_PercentChange_WIDE.csv")
#     with open(csv_wide, "w", encoding="utf-8") as f:
#         f.write("Region,9R,FA\n")
#         for r in regions_order:
#             v9 = last_val(per_group, 'rampp9R', r)
#             vFA = last_val(per_group, 'ramppFA', r)
#             f.write(f"{r},{v9:.6f},{vFA:.6f}\n")

#     print(f"CSV (long): {csv_long}")
#     print(f"CSV (wide): {csv_wide}")


#     # --- Console summary ---
#     def last(x):
#         return x[-1] if x else np.nan

#     print("\n=== SUMMARY (Last step %) ===")
#     print("\nOriginal mapping:")
#     for lbl in [REGION_NAME_MAPPING[r] for r in BRAIN_REGIONS]:
#         v9 = last(per_region['rampp9R'][lbl]); vFA = last(per_region['ramppFA'][lbl])
#         print(f"  {lbl:>3}: 9R={v9:.2f}%   FA={vFA:.2f}%   Diff={v9 - vFA:.2f}%")

#     print("\nGrouped mapping:")
#     for g in GROUPED_REGIONS.keys():
#         v9 = last(per_group['rampp9R'][g]); vFA = last(per_group['ramppFA'][g])
#         print(f"  {g:>3}: 9R={v9:.2f}%   FA={vFA:.2f}%   Diff={v9 - vFA:.2f}%")

#     print(f"\nAll plots saved in: {OUTPUT_DIRECTORY_PATH}")

# if __name__ == "__main__":
#     try:
#         main()
#     except Exception as e:
#         print("ERROR:", e)
#         sys.exit(1)




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate volume-change plots for BOTH mappings:
1) Original REGION_NAME_MAPPING (per-region)
2) Grouped mapping you requested (Am, BG, BS, C, Hi, CC, CR, M, CB)

PLUS: Time evolution plots for each region showing:
  - Absolute volume
  - Volume fraction (region volume / total brain volume)
  - Percentage change in volume fraction

Usage:
    python volume_change_plots_both_mappings.py
"""

import os
import sys
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use('Agg')  # headless
import matplotlib.pyplot as plt

# -------------------- Config --------------------
INPUT_DIRECTORY_PATH = './output_atrophy/post_processing/rampp_final/brainrampp/'
OUTPUT_DIRECTORY_PATH = './output_atrophy/post_processing/rampp_final/plots/volume_evolution_per_region/'

BASE_FILENAMES = ['rampp9R', 'ramppFA']  # models
TIMESTEPS = list(range(17))  # 0..16

# Original mapping setup
BRAIN_REGIONS = [1011, 7, 2, 174, 173, 175, 1030, 17, 1024, 1028, 11, 251, 13, 12, 1035, 10, 18, 4]
REGION_NAME_MAPPING = {
    1024: "MC", 1030: "TL", 1028: "FC", 2: "CR", 1035: "CI", 12: "Pu", 13: "Pa", 10: "Th",
    251: "CC", 11: "NC", 18: "Am", 17: "Hi", 173: "M", 174: "P", 1011: "VC", 7: "CB", 175: "Me", 4: "V"
}

GROUPED_REGIONS = {
    'Corpus Callosum': [251],               # CC = CC
    'Corona Radiata': [2],                 # CR = CR
    'Hippocampus': [17],                # Hi = Hi
    'Cerebellum': [7],                  # CB = CB
    'Midbrain':  [173, 10],           # M = M, Th
    'Brain Stem': [174, 175],          # BS = P, Me
    'Basal Ganglia': [11, 13, 12],        # BG = NC, Pa, Pu   
    'Amygdala': [18],                # Am = Am
    'Cortex':  [1030, 1035, 1028, 1011, 1024],  # C = TL, CI, FC, VC, MC
    'Ventricles': [4],                  # V = V
    "CSF": [24]
}

# Colors (Oxford blue & OI yellow)
COLOR_9R = np.array([0, 114, 178]) / 255.0
COLOR_FA = np.array([240, 228, 66]) / 255.0

os.makedirs(OUTPUT_DIRECTORY_PATH, exist_ok=True)

# -------------------- Helpers --------------------

def compute_hex_volume(coords: np.ndarray) -> float:
    """Approximate hexahedral cell volume from 8 corner points."""
    if len(coords) == 8:
        v0 = coords[0]
        return abs(np.dot(np.cross(coords[1] - v0, coords[2] - v0), coords[4] - v0)) / 6.0
    return 0.0

def read_material_array(ugrid) -> np.ndarray:
    """Return the material id array (point data) as numpy array; raise if missing."""
    material_ids_array = None
    pd = ugrid.GetPointData()
    for idx in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(idx)
        if arr and arr.GetName() and "material" in arr.GetName().lower():
            material_ids_array = arr
            break
    if material_ids_array is None:
        raise RuntimeError("Material IDs array not found in VTK file.")
    return vtk_to_numpy(material_ids_array)

def load_ugrid(vtk_path: str):
    """Read a legacy VTK unstructured grid file."""
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_path)
    reader.Update()
    return reader.GetOutput()

def accumulate_volumes(ugrid, material_ids_np, region_ids):
    """Compute total volume per material id in region_ids and total/24 bookkeeping."""
    volumes = {rid: 0.0 for rid in region_ids}
    volumes.update({'entire_volume': 0.0, 24: 0.0})
    n_cells = ugrid.GetNumberOfCells()
    for cid in range(n_cells):
        cell = ugrid.GetCell(cid)
        pid0 = cell.GetPointIds().GetId(0)
        mat_id = int(material_ids_np[pid0])
        pts = cell.GetPoints()
        coords = np.array([pts.GetPoint(i) for i in range(pts.GetNumberOfPoints())], dtype=float)
        v = compute_hex_volume(coords)
        volumes['entire_volume'] += v
        if mat_id in volumes:
            volumes[mat_id] += v
        # 24 bookkeeping
        if 24 in volumes and mat_id == 24:
            volumes[24] += v
    return volumes

def summarize_last_step_changes(all_data, label_map, outfile_prefix, title_suffix=""):
    """Make grouped bar plots (unsorted and sorted) for last-step percentage changes."""
    # Extract last-step changes & labels
    labels = []
    changes_9R = []
    changes_FA = []
    for key in label_map:  # keep the natural order of label_map
        labels.append(key)
        changes_9R.append(all_data['rampp9R'][key][-1])
        changes_FA.append(all_data['ramppFA'][key][-1])

    # Helper to label bars
    def add_value_labels(ax, bars, values):
        for bar, value in zip(bars, values):
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., 
                    h + (0.5 if h >= 0 else -1.0),
                    f'{value:.1f}%',
                    ha='center', va='bottom' if h >= 0 else 'top', fontsize=9)

    # --- Unsorted ---
    x = np.arange(len(labels))
    bw = 0.35
    fig, ax = plt.subplots(figsize=(16, 10))
    b1 = ax.bar(x - bw/2, changes_9R, bw,
            color=COLOR_9R, edgecolor='black',
            label='9R region based μ', alpha=0.85, linewidth=0.8)

    b2 = ax.bar(x + bw/2, changes_FA, bw,
            color=COLOR_FA, edgecolor='black',
            label='FA voxel based μ', alpha=0.85, linewidth=0.8)
    add_value_labels(ax, b1, changes_9R); add_value_labels(ax, b2, changes_FA)
    ax.set_ylabel("Percentage Change in Volume (%)", fontsize=20)
    ax.set_xticks(x, labels, rotation=45, ha='right', fontsize=14)
    plt.yticks(fontsize=14)
    ax.legend(fontsize=18, loc='lower right')
    fig.tight_layout()
    out_unsorted = os.path.join(OUTPUT_DIRECTORY_PATH, f"{outfile_prefix}_Last_Step_Change.png")
    fig.savefig(out_unsorted, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # --- Sorted by average of 9R and FA ---
    combo = list(zip(labels, changes_9R, changes_FA))
    combo.sort(key=lambda t: (t[1] + t[2]) / 2.0)
    slabels, s9, sFA = zip(*combo) if combo else ([], [], [])
    x = np.arange(len(slabels))
    fig, ax = plt.subplots(figsize=(16, 10))
    b1 = ax.bar(x - bw/2, s9, bw,
                color=COLOR_9R, edgecolor='black',
                label='9R region based μ', alpha=0.85, linewidth=0.8)

    b2 = ax.bar(x + bw/2, sFA, bw,
                color=COLOR_FA, edgecolor='black',
                label='FA voxel based μ', alpha=0.85, linewidth=0.8)
    add_value_labels(ax, b1, s9); add_value_labels(ax, b2, sFA)
    ax.set_xlabel("Brain Regions" if title_suffix == "Original" else "Brain Regions", fontsize=20)
    ax.set_ylabel("Percentage Change at Last Step (%)", fontsize=20)
    ax.set_xticks(x, slabels, rotation=45, ha='right', fontsize=14)
    plt.yticks(fontsize=14)
    ax.legend(fontsize=18, loc='lower right')
    fig.tight_layout()
    out_sorted = os.path.join(OUTPUT_DIRECTORY_PATH, f"{outfile_prefix}_Last_Step_Change_Sorted.png")
    fig.savefig(out_sorted, dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved: {out_unsorted}")
    print(f"Saved: {out_sorted}")

def plot_time_evolution(time_data, region_mapping, outfile_prefix, title_suffix=""):
    """
    Create time evolution plots for each region showing:
    1. Absolute volume over time
    2. Volume fraction over time
    3. Percentage change in volume fraction
    
    time_data structure: {base: {region: {'volume': [], 'brain_vol': [], 'timesteps': []}}}
    """
    print(f"\nGenerating time evolution plots for {title_suffix} mapping...")
    
    for region in region_mapping:
        fig, axes = plt.subplots(3, 1, figsize=(12, 14))
        fig.suptitle(f'Time Evolution - {region}', fontsize=16, fontweight='bold')
        
        for base in BASE_FILENAMES:
            data = time_data[base][region]
            timesteps = data['timesteps']
            volumes = data['volume']
            brain_vols = data['brain_vol']
            
            # Skip if no data
            if not timesteps:
                continue
            
            # Calculate volume fractions
            vol_fractions = [v / bv if bv > 0 else 0 for v, bv in zip(volumes, brain_vols)]
            
            # Calculate percentage changes in volume fraction (relative to first timestep)
            if vol_fractions:
                initial_frac = vol_fractions[0]
                pct_changes = [100 * (vf - initial_frac) / initial_frac if initial_frac > 0 else 0 
                              for vf in vol_fractions]
            else:
                pct_changes = []
            
            color = COLOR_9R if base == 'rampp9R' else COLOR_FA
            label = '9R region based μ' if base == 'rampp9R' else 'FA voxel based μ'
            
            # Plot 1: Absolute Volume
            axes[0].plot(timesteps, volumes, marker='o', linewidth=2, 
                        color=color, label=label, markersize=5)
            
            # Plot 2: Volume Fraction
            axes[1].plot(timesteps, vol_fractions, marker='o', linewidth=2,
                        color=color, label=label, markersize=5)
            
            # Plot 3: Percentage Change in Volume Fraction
            axes[2].plot(timesteps, pct_changes, marker='o', linewidth=2,
                        color=color, label=label, markersize=5)
        
        # Configure Plot 1: Absolute Volume
        axes[0].set_ylabel('Volume (mm³)', fontsize=12, fontweight='bold')
        axes[0].grid(True, alpha=0.3)
        axes[0].legend(fontsize=10, loc='best')
        axes[0].set_title('Absolute Volume', fontsize=12)
        
        # Configure Plot 2: Volume Fraction
        axes[1].set_ylabel('Volume Fraction', fontsize=12, fontweight='bold')
        axes[1].grid(True, alpha=0.3)
        axes[1].legend(fontsize=10, loc='best')
        axes[1].set_title('Volume Fraction (Region/Brain)', fontsize=12)
        
        # Configure Plot 3: Percentage Change
        axes[2].set_xlabel('Timestep', fontsize=12, fontweight='bold')
        axes[2].set_ylabel('Change in Volume Fraction (%)', fontsize=12, fontweight='bold')
        axes[2].grid(True, alpha=0.3)
        axes[2].legend(fontsize=10, loc='best')
        axes[2].set_title('Percentage Change in Volume Fraction', fontsize=12)
        axes[2].axhline(y=0, color='k', linestyle='--', alpha=0.3)
        
        fig.tight_layout()
        
        # Save figure
        safe_region_name = region.replace(' ', '_').replace('/', '_')
        outfile = os.path.join(OUTPUT_DIRECTORY_PATH, 
                              f"{outfile_prefix}_TimeEvolution_{safe_region_name}.png")
        fig.savefig(outfile, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {outfile}")

# -------------------- Main processing --------------------
def main():
    # Structures to hold % changes per timestep for ORIGINAL region labels (by display name)
    per_region = {base: {REGION_NAME_MAPPING[r]: [] for r in BRAIN_REGIONS} for base in BASE_FILENAMES}
    # Structures to hold % changes per timestep for GROUPED labels
    per_group = {base: {g: [] for g in GROUPED_REGIONS.keys()} for base in BASE_FILENAMES}
    
    # NEW: Structures to hold time evolution data
    time_data_original = {
        base: {
            REGION_NAME_MAPPING[r]: {'volume': [], 'brain_vol': [], 'timesteps': []} 
            for r in BRAIN_REGIONS
        } for base in BASE_FILENAMES
    }
    
    time_data_grouped = {
        base: {
            g: {'volume': [], 'brain_vol': [], 'timesteps': []} 
            for g in GROUPED_REGIONS.keys()
        } for base in BASE_FILENAMES
    }

    for base in BASE_FILENAMES:
        print(f"\nProcessing model: {base}")
        # initial absolute volumes for each ORIGINAL region / GROUP for step 0
        init_region_vols = None
        init_group_vols = None

        for step in TIMESTEPS:
            vtk_path = os.path.join(INPUT_DIRECTORY_PATH, f"{base}_{step}.vtk")
            if not os.path.exists(vtk_path):
                print(f"  Missing file: {vtk_path} (skipping this step)")
                continue

            print(f"  Step {step}: {vtk_path}")
            ugrid = load_ugrid(vtk_path)
            mat_np = read_material_array(ugrid)

            # --- ORIGINAL region volumes at this step ---
            vols_region = accumulate_volumes(ugrid, mat_np, BRAIN_REGIONS)
            # Brain volume excluding 24
            brain_volume = vols_region['entire_volume'] - vols_region.get(24, 0.0)

            # Track absolute per-region
            abs_region = {REGION_NAME_MAPPING[r]: vols_region.get(r, 0.0) for r in BRAIN_REGIONS}

            # --- GROUPED volumes at this step ---
            abs_group = {}
            for gname, members in GROUPED_REGIONS.items():
                abs_group[gname] = float(sum(vols_region.get(mid, 0.0) for mid in members))

            # Initialize baselines at first encountered step
            if init_region_vols is None:
                init_region_vols = abs_region.copy()
            if init_group_vols is None:
                init_group_vols = abs_group.copy()

            # Store time evolution data for ORIGINAL regions
            for lbl in abs_region:
                time_data_original[base][lbl]['volume'].append(abs_region[lbl])
                time_data_original[base][lbl]['brain_vol'].append(brain_volume)
                time_data_original[base][lbl]['timesteps'].append(step)
            
            # Store time evolution data for GROUPED regions
            for g in abs_group:
                time_data_grouped[base][g]['volume'].append(abs_group[g])
                time_data_grouped[base][g]['brain_vol'].append(brain_volume)
                time_data_grouped[base][g]['timesteps'].append(step)

            # Percentage changes relative to step-0
            for lbl in abs_region:
                base_v = init_region_vols[lbl]
                delta_pct = 0.0 if base_v == 0 else 100.0 * (abs_region[lbl] - base_v) / base_v
                per_region[base][lbl].append(delta_pct)

            for g in abs_group:
                base_v = init_group_vols[g]
                delta_pct = 0.0 if base_v == 0 else 100.0 * (abs_group[g] - base_v) / base_v
                per_group[base][g].append(delta_pct)

    # --- Plots for ORIGINAL mapping ---
    print("\nGenerating plots for ORIGINAL mapping (last-step %) ...")
    summarize_last_step_changes(per_region, 
                                label_map=[REGION_NAME_MAPPING[r] for r in BRAIN_REGIONS],
                                outfile_prefix="OriginalMapping",
                                title_suffix="Original")

    # --- Plots for GROUPED mapping ---
    print("\nGenerating plots for GROUPED mapping (last-step %) ...")
    summarize_last_step_changes(per_group, 
                                label_map=list(GROUPED_REGIONS.keys()),
                                outfile_prefix="GroupedMapping",
                                title_suffix="Grouped")
    
    # --- NEW: Time evolution plots for ORIGINAL mapping ---
    plot_time_evolution(time_data_original, 
                       [REGION_NAME_MAPPING[r] for r in BRAIN_REGIONS],
                       "OriginalMapping",
                       "Original")
    
    # --- NEW: Time evolution plots for GROUPED mapping ---
    plot_time_evolution(time_data_grouped,
                       list(GROUPED_REGIONS.keys()),
                       "GroupedMapping",
                       "Grouped")

    # --- Export CSV tables for LaTeX/pgfplots (9 grouped regions, last-step %) ---
    regions_order = list(GROUPED_REGIONS.keys())

    def last_val(dct, base, key):
        vals = dct[base][key]
        return vals[-1] if vals else np.nan

    # Long format: Region,Model,PercentChange
    csv_long = os.path.join(OUTPUT_DIRECTORY_PATH, "GroupedMapping_LastStep_PercentChange_LONG.csv")
    with open(csv_long, "w", encoding="utf-8") as f:
        f.write("Region,Model,PercentChange\n")
        for r in regions_order:
            f.write(f"{r},9R,{last_val(per_group, 'rampp9R', r):.6f}\n")
            f.write(f"{r},FA,{last_val(per_group, 'ramppFA', r):.6f}\n")

    # Wide format: Region,9R,FA
    csv_wide = os.path.join(OUTPUT_DIRECTORY_PATH, "GroupedMapping_LastStep_PercentChange_WIDE.csv")
    with open(csv_wide, "w", encoding="utf-8") as f:
        f.write("Region,9R,FA\n")
        for r in regions_order:
            v9 = last_val(per_group, 'rampp9R', r)
            vFA = last_val(per_group, 'ramppFA', r)
            f.write(f"{r},{v9:.6f},{vFA:.6f}\n")

    print(f"\nCSV (long): {csv_long}")
    print(f"CSV (wide): {csv_wide}")

    # --- Console summary ---
    def last(x):
        return x[-1] if x else np.nan

    print("\n=== SUMMARY (Last step %) ===")
    print("\nOriginal mapping:")
    for lbl in [REGION_NAME_MAPPING[r] for r in BRAIN_REGIONS]:
        v9 = last(per_region['rampp9R'][lbl]); vFA = last(per_region['ramppFA'][lbl])
        print(f"  {lbl:>3}: 9R={v9:.2f}%   FA={vFA:.2f}%   Diff={v9 - vFA:.2f}%")

    print("\nGrouped mapping:")
    for g in GROUPED_REGIONS.keys():
        v9 = last(per_group['rampp9R'][g]); vFA = last(per_group['ramppFA'][g])
        print(f"  {g:>3}: 9R={v9:.2f}%   FA={vFA:.2f}%   Diff={v9 - vFA:.2f}%")

    print(f"\nAll plots saved in: {OUTPUT_DIRECTORY_PATH}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR:", e)
        sys.exit(1)