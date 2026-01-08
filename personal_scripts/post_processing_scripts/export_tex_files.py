import numpy as np
import pandas as pd
import vtk
import os
import seaborn as sns
import matplotlib.pyplot as plt
from os.path import exists
from matplotlib.lines import Line2D
import scipy.stats as sp

# --- CONFIGURATION ---
base_path = "./output_atrophy/post_processing/rampp_final/"
fileNames = [
    "rampp_atrophy_9R_last_timestep",
    "rampp_atrophy_FA_last_timestep"
]
output_dir = os.path.join(base_path, "boxplots_stats")
os.makedirs(output_dir, exist_ok=True)

measure = [
    'displacement', 'max_principal_stretch', 'med_principal_stretch',
    'min_principal_stretch', 'max_shear', 'hyrdostatic_stress', 'von_mises'
]

# measure = [
#     'displacement'
# ]

# measure = [
#     'hyrdostatic_stress', 'von_mises'
# ]

metric_labels = {
    'max_principal_stretch': 'First principal stretch [-]',
    'med_principal_stretch': 'Second principal stretch [-]',
    'min_principal_stretch': 'Third principal stretch [-]',
    'max_shear': 'Maximum shear [-]',
    'hyrdostatic_stress': 'Hydrostatic stresses [kPa]',
    'von_mises': 'von Mises stresses [kPa]',
    'displacement': 'Displacement [mm]'
}

# Keep this exact order (matches your LaTeX xticks)
regions = {
    'Corpus Callosum': [251],
    'Corona Radiata': [2],
    'Hippocampus': [17],
    'Cerebellum': [7],
    'Midbrain':  [173, 10, 30000],
    'Brain Stem': [174, 175],
    'Basal Ganglia': [11, 13, 12, 18000],
    'Amygdala': [18],
    'Cortex':  [1030, 1035, 1028, 1011, 1024]}

# regions = {
#     'Corpus Callosum': [251,252,253,254,255],
#     'Corona Radiata': [2,41,77,78,79,1026,1002,1023,1010,2026,2002,2023,2010],
#     'Hippocampus': [17,53],
#     'Cerebellum': [7,46,8,47, 60000],
#     'Midbrain':  [173,10,49,28,60, 24000],

#     'Brain Stem': [174, 175, 12000],
#     'Basal Ganglia': [11, 12, 13, 51,26,58,52,50],
#     'Amygdala': [18, 54, 9000],
#     'Cortex':  [1030,1015,1009,1007,1034,1006,1016,
#                        2030,2015,2009,2007,2034,2006,2016,
#                        1035,2035,1028,1027,1003,1018,1020,
#                        1019,1012,1014,1017,2028,2027,2003,
#                        2018,2020,2019,2012,2014,2017,1011,
#                        1013,1005,1021,2011,2013,2005,2021,
#                        1024,1029,1008,1031,1022,1025,2024,
#                        2029,2008,2031,2022,2025, 6000],
# }

region_order = [r for r in regions.keys() if r != 'Ventricles']
# region_id_to_label = {v[0]: k for k, v in regions.items()}

region_id_to_label = {}
overlaps = {}
for label, ids in regions.items():
    for rid in ids:
        if rid in region_id_to_label and region_id_to_label[rid] != label:
            overlaps.setdefault(rid, set()).update({region_id_to_label[rid], label})
        region_id_to_label[rid] = label


# Sanity checks (optional, but very useful)
present_ids = set()  # fill after you read a file; shown here as a function for reuse
def validate_mapping(material_ids_in_mesh):
    present = set(material_ids_in_mesh)
    defined = set(region_id_to_label.keys())
    unused_defined = sorted(defined - present)
    undefined_present = sorted(present - defined)

    if overlaps:
        print(f"[WARN] Some material IDs are assigned to multiple regions: "
              + ", ".join(f"{rid}:{'/'.join(sorted(labs))}" for rid, labs in overlaps.items()))
    if unused_defined:
        print(f"[INFO] Region map contains IDs not in mesh (ignored): {unused_defined[:20]}"
              + (" ..." if len(unused_defined) > 20 else ""))
    if undefined_present:
        print(f"[WARN] Mesh has material IDs with no region mapping (dropped): {undefined_present[:20]}"
              + (" ..." if len(undefined_present) > 20 else ""))

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
    pd = data.GetPointData()
    cd = data.GetCellData()

    # --- find material IDs (point pref, else derive from cell data) ---
    material_ids = None
    # 1) try point data
    for a in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(a)
        if not arr or not arr.GetName():
            continue
        if "material" in arr.GetName().lower():
            material_ids = [int(arr.GetComponent(i, 0)) for i in range(num_points)]
            break


    # 2) else try cell data and project to points via majority of incident cells
    if material_ids is None:
        cell_mat_arr = None
        for a in range(cd.GetNumberOfArrays()):
            arr = cd.GetArray(a)
            if arr and arr.GetName() and "material" in arr.GetName().lower():
                cell_mat_arr = arr
                break
        if cell_mat_arr is not None:
            idlist = vtk.vtkIdList()
            material_ids = []
            for i in range(num_points):
                data.GetPointCells(i, idlist)
                if idlist.GetNumberOfIds() == 0:
                    material_ids.append(-1)
                    continue
                vals = [int(cell_mat_arr.GetComponent(idlist.GetId(k), 0))
                        for k in range(idlist.GetNumberOfIds())]
                # majority vote
                material_ids.append(max(set(vals), key=vals.count))
        else:
            print(f"No material IDs found in {fullFilename}. Skipping.")
            return None

    validate_mapping(material_ids)
    # --- collect arrays; compute displacement magnitude if needed ---
    filenameData = {}  # name -> list of per-point scalar values

    for a in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(a)
        if not arr or not arr.GetName():
            continue
        name = arr.GetName()
        lname = name.lower()
        if "material" in lname:
            continue

        ncomp = arr.GetNumberOfComponents()

        # displacement vector -> magnitude
        if ncomp == 3 and ("disp" in lname or "displacement" in lname):
            vals = []
            for i in range(num_points):
                t = arr.GetTuple(i)  # (ux, uy, uz)
                mag = float((t[0]**2 + t[1]**2 + t[2]**2) ** 0.5)
                vals.append(mag)
            filenameData["displacement"] = vals  # canonical key
        elif ncomp == 1:
            vals = [float(arr.GetComponent(i, 0)) for i in range(num_points)]
            filenameData[name] = vals
        else:
            continue

    if len(filenameData) == 0:
        print(f"[WARN] No usable scalar arrays in {fullFilename}.")
        return None

    # --- group by region/material id ---
    region_data = {}  # metric -> {region_id: [values]}
    for key, array in filenameData.items():
        region_values = {}
        for i, mat_id in enumerate(material_ids):
            region_values.setdefault(mat_id, []).append(array[i])
        region_data[key if key in measure else key.lower()] = region_values

    # displacement must be >= 0
    if "displacement" in region_data:
        for rid, vals in region_data["displacement"].items():
            region_data["displacement"][rid] = [v if v >= 0.0 else 0.0 for v in vals]

    return region_data

def remove_outliers_percentile(values, lower=5, upper=95):
    """Remove outliers from a list of values using percentile method."""
    if not values:
        return values
    arr = np.array(values)
    
    # For very small ranges, be more conservative with outlier removal
    value_range = np.max(arr) - np.min(arr)
    if value_range < 1e-6:  # Very small range, don't remove outliers
        return arr.tolist()
    
    low = np.percentile(arr, lower)
    high = np.percentile(arr, upper)
    
    # Ensure we don't remove too much data for small ranges
    filtered = arr[(arr >= low) & (arr <= high)]
    if len(filtered) < len(arr) * 0.5:  # If we'd remove more than 50% of data, be less aggressive
        low = np.percentile(arr, 1)  # Use 1st and 99th percentile instead
        high = np.percentile(arr, 99)
        filtered = arr[(arr >= low) & (arr <= high)]
    
    return filtered.tolist()

def compute_wilcoxon_cohens_d(arr1, arr2):
    arr1 = np.array(arr1)
    arr2 = np.array(arr2)
    min_len = min(len(arr1), len(arr2))
    if min_len == 0:
        return np.nan, np.nan, np.nan
    if len(arr1) != len(arr2):
        arr1 = arr1[:min_len]
        arr2 = arr2[:min_len]
    differences = arr1 - arr2
    mean_diff = np.mean(differences)

    if mean_diff > 0:
        alternative = 'greater'
    elif mean_diff < 0:
        alternative = 'less'
    else:
        alternative = 'two-sided'

    try:
        non_zero_diffs = differences[differences != 0]
        if len(non_zero_diffs) == 0:
            wilcoxon_stat, wilcoxon_p = 0.0, 1.0
        else:
            wilcoxon_stat, wilcoxon_p = sp.wilcoxon(non_zero_diffs, alternative=alternative)
    except ValueError:
        wilcoxon_stat, wilcoxon_p = np.nan, np.nan

    std_diff = np.std(differences, ddof=1)
    cohens_d = (mean_diff / std_diff) if std_diff > 0 else np.nan
    return wilcoxon_stat, wilcoxon_p, cohens_d

# --- DATA LOADING AND ORGANIZATION ---
all_data = []
model_names = []
data_by_file = []
for f in fileNames:
    full_path = os.path.join(base_path, f + ".vtk")
    region_data = load_vtk_data(full_path)
    if region_data is None:
        continue
    data_by_file.append(region_data)
    model_name = f.split('_')[-3]
    model_names.append(model_name)

    for metric in measure:
        if metric not in region_data:
            continue
        for region_id, values in region_data[metric].items():
            # Skip ventricles & any unknown region
            if region_id in [4, 24]:
                continue
            region_label = region_id_to_label.get(region_id)
            if region_label is None:
                continue

            # --- unit conversion for stresses ---
            if metric in ["hyrdostatic_stress", "von_mises"]:
                values = [v * 1000.0 for v in values]  # MPa → kPa

            for v in values:
                all_data.append({
                    'Region': region_label,
                    'Value': v,
                    'Atrophy model': model_name,
                    'Metric': metric
                })

def compute_boxplot_summary(values):
    """Return dict compatible with PGFPlots 'boxplot prepared={...}'."""
    arr = np.asarray(values)
    if arr.size == 0:
        return None
    q5  = float(np.percentile(arr, 5))
    q25 = float(np.percentile(arr, 25))
    q50 = float(np.percentile(arr, 50))
    q75 = float(np.percentile(arr, 75))
    q95 = float(np.percentile(arr, 95))
    mean = float(np.mean(arr))
    return {
        "lower_whisker": q5,
        "lower_quartile": q25,
        "median": q50,
        "upper_quartile": q75,
        "upper_whisker": q95,
        "mean": mean,
        "max": float(np.max(arr)),
        "min": float(np.min(arr))
    }

def make_xtick_labels_with_diamonds():
    """Build xticklabels with \\diamondsymbol colors per region."""
    # region_colors = {
    #     'Corpus Callosum': 'OI-Vermilion',
    #     'Corona Radiata': 'OI-White',
    #     'Hippocampus': 'OI-Orange',
    #     'Cerebellum': 'OI-Green',
    #     'Midbrain': 'OI-Yellow',
    #     'Brain Stem': 'OI-Blue',
    #     'Basal Ganglia': 'OI-Purple',
    #     'Amygdala': 'lime',
    #     'Cortex': 'OI-SkyBlue',
    # }

    region_colors = {
        'Corpus Callosum': 'OI-Red',
        'Corona Radiata': 'OI-Gray',
        'Hippocampus': 'OI-Purple',
        'Cerebellum': 'OI-LightGreen',
        'Midbrain': 'OI-DarkBlue',
        'Brain Stem': 'OI-Darkred',
        'Basal Ganglia': 'OI-Yellow',
        'Amygdala': 'OI-Brown',
        'Cortex': 'OI-DarkGreen',
    }
    parts = []
    for r in region_order:
        color = region_colors.get(r, 'OI-Black')
        parts.append(f"\\diamondsymbol{{{color}}}~{r}")
    return ",\n          ".join(parts)

# def tex_sanitize_float(x):
#     """Guard against nan/inf in TeX numbers."""
#     if x is None or not np.isfinite(x):
#         return "nan"
#     return f"{x:.3f}"

def tex_sanitize_float(x, digits=8):
    if x is None or not np.isfinite(x):
        return "nan"
    return f"{float(x):.{digits}f}"

def effect_color_and_text(d):
    """Return (tikz_color_name, legend_text) for a Cohen's d value."""
    if not np.isfinite(d):
        return ("OI-Black", "$d=\\mathrm{NA}$")
    ad = abs(d)
    if ad < 0.499:
        return ("OI-Green", f"$d={d:.2f}$")   # Small
    elif ad < 0.799:
        return ("OI-Orange", f"$d={d:.2f}$")  # Medium
    else:
        return ("OI-Red", f"$d={d:.2f}$")     # Large

def write_tex_figure_for_metric(df_all: pd.DataFrame, metric: str):
    """Emit a LaTeX figure with Cohen's d labels and outlier-filtered y-range."""
    df = df_all[df_all['Metric'] == metric].copy()
    if df.empty:
        print(f"[WARN] Skipping TeX for {metric}: no data.")
        return

    # Map models to Sim1/Sim2 by input order
    if len(model_names) < 2:
        print(f"[WARN] Skipping TeX for {metric}: need two models.")
        return
    model_to_sim = {model_names[0]: 'Sim1', model_names[1]: 'Sim2'}
    df['sim'] = df['Atrophy model'].map(model_to_sim)
    df = df[df['Region'].isin(region_order)]

    # Remove outliers from data for y-range calculation
    filtered_values = []
    summaries = {}   # (region, sim) -> dict with quartiles/whiskers/mean/max
    region_ymax = {} # region -> max across both sims (for label placement)
    
    for r in region_order:
        ymax_r = None
        for sim in ('Sim1', 'Sim2'):
            vals = df[(df['Region'] == r) & (df['sim'] == sim)]['Value'].to_numpy()
            # Remove outliers for this region/sim combination
            vals_filtered = remove_outliers_percentile(vals.tolist(), lower=5, upper=95)
            filtered_values.extend(vals_filtered)
            
            summ = compute_boxplot_summary(vals_filtered)  # Use filtered data for summary
            summaries[(r, sim)] = summ
            if summ:
                ymax_r = summ['max'] if ymax_r is None else max(ymax_r, summ['max'])
        region_ymax[r] = ymax_r if ymax_r is not None else 0.0

    # Compute Cohen's d per region from original (non-filtered) data for statistical validity
    def compute_local_d(region_name):
        sim1 = df[(df['Region'] == region_name) & (df['sim'] == 'Sim1')]['Value'].to_numpy()
        sim2 = df[(df['Region'] == region_name) & (df['sim'] == 'Sim2')]['Value'].to_numpy()
        if sim1.size == 0 or sim2.size == 0:
            return np.nan
        n = min(sim1.size, sim2.size)
        sim1 = sim1[:n]
        sim2 = sim2[:n]
        diffs = sim1 - sim2
        sd = np.std(diffs, ddof=1)
        if sd == 0 or not np.isfinite(sd):
            return 0.0
        d = float(np.mean(diffs) / sd)
        return d

    # Calculate y-limits based on filtered data
    if filtered_values:
        ymax_filtered = max(filtered_values)
        ymin_filtered = min(filtered_values)
        y_range = ymax_filtered - ymin_filtered
        
        # Handle cases where range is very small or zero
        if y_range < 1e-10:  # Very small range
            y_center = (ymax_filtered + ymin_filtered) / 2
            y_range = max(abs(ymax_filtered), abs(ymin_filtered)) * 0.1  # 10% of max absolute value
            if y_range < 1e-10:  # Still too small, use default
                y_range = 1e-6
            ymin_default = y_center - y_range
            ymax_default = y_center + y_range
        else:
            # Normal case with reasonable range
            y_headroom = max(0.15 * y_range, abs(ymax_filtered) * 0.05)  # At least 15% headroom or 5% of max value
            ymax_default = ymax_filtered + y_headroom
            ymin_default = ymin_filtered - 0.02 * y_range
    else:
        ymax_default = 1.0
        ymin_default = 0.0

    # Axis limits and tick positions
    base_positions = [i*3 for i in range(len(region_order))]  # 0,3,6,...
    xticks = [b+1 for b in base_positions]                    # 1,4,7,...
    xticklabels = make_xtick_labels_with_diamonds()

    ylabel = metric_labels.get(metric, metric)

    tex_path = os.path.join(output_dir, f"{metric}.tex")
    with open(tex_path, 'w') as f:
        f.write("% Auto-generated by Python: boxplots with Cohen's d labels (outlier-filtered y-range)\n")
        f.write("%%%%%%%%%%%%%%%%%%%%%%%%\n")
        f.write(f"%Diagram {metric}\n")
        f.write("%%%%%%%%%%%%%%%%%%%%%%\n\n\n")

        f.write("\\begin{figure*}[p]\n")
        f.write("  \\centering\n")
        f.write("  \\begin{tikzpicture}\n")
        f.write("    \\begin{axis}[\n")
        f.write("      width=\\textwidth,\n")
        f.write("      height=0.4\\textheight,\n")
        f.write(f"      ymin={tex_sanitize_float(ymin_default)}, ymax={tex_sanitize_float(ymax_default)},\n")
        f.write("      xtick={%s},\n" % ",".join(str(x) for x in xticks))
        f.write("      xticklabels={\n          %s\n      },\n" % xticklabels)
        f.write("      legend to name=fulllegend,\n")
        f.write("      legend style={\n")
        f.write("        draw=black,\n")
        f.write("        fill=white,\n")
        f.write("        legend cell align=left,\n")
        f.write("        legend columns=3,\n")
        f.write("        column sep=0.4cm\n")
        f.write("      },\n")
        f.write("      x tick label style={rotate=45,anchor=east,font=\\small},\n")
        f.write("      xmin=0.5, xmax=%s,\n" % tex_sanitize_float(base_positions[-1] + 2.5))
        f.write(f"      ylabel={{{{ {ylabel} }}}},\n")
        f.write("      xlabel={Brain regions},\n")
        f.write("      enlarge x limits=0.15,\n")
        f.write("      boxplot/draw direction=y,\n")
        f.write("      boxplot/remove outliers,\n")
        f.write("      boxplot/whisker line style={solid},\n")
        f.write("      boxplot/median style={solid},\n")
        f.write("      every box/.style={solid},\n")
        f.write("      every average/.style={only marks, mark=triangle, mark options={fill=white}, solid},\n")
        f.write("    ]\n\n")

        # Legend entries
        # f.write("      % manual legend entries\n")
        # f.write("      \\addlegendimage{area legend, draw=OI-Black, fill=OI-Blue, opacity=0.8}\n")
        # f.write("      \\addlegendentry{boxplot region based $\\mu$}\n\n")
        # f.write("      \\addlegendimage{area legend, draw=OI-Black, fill=OI-Yellow, opacity=1}\n")
        # f.write("      \\addlegendentry{boxplot voxel based $\\mu$}\n\n")
        # f.write("      \\addlegendimage{only marks, mark=o, mark size=3pt, draw=OI-Black, fill=white}\n")
        # f.write("      \\addlegendentry{arithmetic mean}\n\n")
        # f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Green}\n")
        # f.write("      \\addlegendentry{Small effect ($d<0.499$)}\n\n")
        # f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Orange}\n")
        # f.write("      \\addlegendentry{Medium effect ($0.5\\le d<0.799$)}\n\n")
        # f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Red}\n")
        # f.write("      \\addlegendentry{Large effect ($d\\ge0.8$)}\n\n")


        # f.write("      % manual legend entries\n")
        # f.write("      \\addlegendimage{area legend, draw=OI-Black, fill=OI-Blue, opacity=0.8}\n")
        # f.write("      \\addlegendentry{boxplot region based $\\mu$}\n\n")
        # f.write("      \\addlegendimage{area legend, draw=OI-Black, fill=OI-Yellow, opacity=1}\n")
        # f.write("      \\addlegendentry{boxplot voxel based $\\mu$}\n\n")
        # f.write("      \\addlegendimage{only marks, mark=o, mark size=3pt, draw=OI-Black, fill=white}\n")
        # f.write("      \\addlegendentry{arithmetic mean}\n\n")
        # f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Green}\n")
        # f.write("      \\addlegendentry{Small effect ($-0.499<d<0.499$)}\n\n")
        # f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Orange}\n")
        # f.write("      \\addlegendentry{Medium effect ($-0.799<d\\le-0.5$ and $0.5\\le d<0.799$)}\n\n")
        # f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Red}\n")
        # f.write("      \\addlegendentry{Large effect ($-0.8\\led$ and $d\\ge0.8$)}\n\n")


        f.write("      % manual legend entries\n")

        f.write("      \\addlegendimage{area legend, draw=OI-Black, fill=OI-Blue, opacity=0.8}\n")
        f.write("      \\addlegendentry{boxplot region based $\\mu$}\n\n")

        f.write("      \\addlegendimage{area legend, draw=OI-Black, fill=OI-Yellow, opacity=1}\n")
        f.write("      \\addlegendentry{boxplot voxel based $\\mu$}\n\n")

        f.write("      \\addlegendimage{only marks, mark=o, mark size=3pt, draw=OI-Black, fill=white}\n")
        f.write("      \\addlegendentry{arithmetic mean}\n\n")

        f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Green}\n")
        f.write("      \\addlegendentry{\\shortstack{Small effect\\\\$-0.499 < d < 0.499$}}\n\n")

        f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Orange}\n")
        f.write("      \\addlegendentry{\\shortstack{Medium effect\\\\$-0.799 < d \\le -0.5$ and $0.5 \\le d < 0.799$}}\n\n")

        f.write("      \\addlegendimage{only marks, mark=square*, mark size=3pt, draw=OI-Black, fill=OI-Red}\n")
        f.write("      \\addlegendentry{\\shortstack{Large effect\\\\$-0.8 < d$ and $d \\ge 0.8$}}\n\n")





        # Emit plots for each region
        for idx, region in enumerate(region_order):
            base = idx*3
            pos_sim1 = base + 0.5
            pos_sim2 = base + 1.5

            for sim, pos, fill in (('Sim1', pos_sim1, 'OI-Blue'), ('Sim2', pos_sim2, 'OI-Yellow')):
                summ = summaries.get((region, sim))
                if summ is None:
                    continue
                uw1 = summaries.get((region, 'Sim1'))
                uw2 = summaries.get((region, 'Sim2'))
                top_whisker = -np.inf
                if uw1 is not None:
                    top_whisker = max(top_whisker, uw1['upper_whisker'])
                if uw2 is not None:
                    top_whisker = max(top_whisker, uw2['upper_whisker'])
                if not np.isfinite(top_whisker):
                    top_whisker = 0.0
                f.write(f"      %--- {region} ({sim}) ---\n")
                f.write("      \\addplot+[\n")
                f.write("        boxplot prepared={\n")
                f.write(f"          lower whisker={tex_sanitize_float(summ['lower_whisker'])},\n")
                f.write(f"          lower quartile={tex_sanitize_float(summ['lower_quartile'])},\n")
                f.write(f"          median={tex_sanitize_float(summ['median'])},\n")
                f.write(f"          upper quartile={tex_sanitize_float(summ['upper_quartile'])},\n")
                f.write(f"          upper whisker={tex_sanitize_float(summ['upper_whisker'])}\n")
                f.write("        },\n")
                f.write(f"        fill={fill}, fill opacity={'0.8' if sim=='Sim1' else '1.0'},\n")
                f.write("        draw=OI-Black, solid,\n")
                f.write(f"        boxplot/draw position={tex_sanitize_float(pos)}\n")
                f.write("      ] coordinates {};\n\n")

                # Mean marker at (x, mean)
                f.write("      \\addplot+[\n")
                f.write("        only marks, mark=o, mark size=3pt,\n")
                f.write("        draw=OI-Black, fill=white,\n")
                f.write(f"        boxplot/draw position={tex_sanitize_float(pos)}\n")
                f.write("      ] coordinates {\n")
                f.write(f"        ({tex_sanitize_float(pos)}, {tex_sanitize_float(summ['mean'])})\n")
                f.write("      };\n\n")

            # Cohen's d label (centered at tick position base+1)
            d_val = compute_local_d(region)
            color_name, d_text = effect_color_and_text(d_val)

            # Position Cohen's d label appropriately above the data
            ylab = region_ymax.get(region, 0.0)
            if not np.isfinite(ylab):
                ylab = 0.0
            
            # Calculate label position based on actual y-range
            y_range_for_label = ymax_default - ymin_default
            y_label_offset = max(0.08 * y_range_for_label, abs(ylab) * 0.1)  # At least 8% of range or 10% of max value
            y_label = ylab + y_label_offset
            x_label = base + 1.0

            # f.write(
            #     f"      \\node[anchor=south, text={color_name}] "
            #     f"at (axis cs:{tex_sanitize_float(x_label)},{tex_sanitize_float(y_label)}) "
            #     f"{{{d_text}}};\n\n"
            # )
            f.write(
                f"      \\node[anchor=south, text={color_name}, yshift=3pt] "
                f"at (axis cs:{tex_sanitize_float(x_label)},{tex_sanitize_float(top_whisker)}) "
                f"{{{d_text}}};\n\n"
            )

        f.write("    \\end{axis}\n\n")
        f.write("    % legend under the axis\n")
        f.write("    \\node[anchor=north]\n")
        f.write("      at (current bounding box.south)\n")
        f.write("      {\\pgfplotslegendfromname{fulllegend}};\n")
        f.write("  \\end{tikzpicture}\n")
        f.write(f"  \\caption{{Boxplots illustrating '{metric}' across nine brain regions. The lower whisker is the 5\\textsuperscript{{th}} percentile and the upper whisker the 95\\textsuperscript{{th}} percentile. The arithmetic mean is marked. Cohen's $d$ is shown above each region.}}\n")
        f.write(f"  \\label{{fig:boxplot_{metric}}}\n")
        f.write("\\end{figure*}\n")

    print(f"[OK] Wrote TeX figure: {tex_path}")

def export_metric_files(df_all: pd.DataFrame, metric: str):
    """Export CSV files for each metric."""
    df = df_all[df_all['Metric'] == metric].copy()
    if df.empty:
        print(f"[WARN] No data for metric: {metric}")
        return

    # Map models to Sim1/Sim2 by input order
    if len(model_names) < 2:
        raise RuntimeError("Need two models to create paired comparisons.")
    model_to_sim = {model_names[0]: 'Sim1', model_names[1]: 'Sim2'}
    df['sim'] = df['Atrophy model'].map(model_to_sim)

    # keep only the 9 regions in the plotting order
    df = df[df['Region'].isin(region_order)]
    df['Region'] = pd.Categorical(df['Region'], categories=region_order, ordered=True)
    df.sort_values(['Region', 'sim'], inplace=True)

    # write metric-specific value column
    value_col = metric
    df_out = df[['Region', 'sim', 'Value']].rename(columns={'Region': 'region', 'Value': value_col})
    data_csv = os.path.join(output_dir, f"{metric}_data.csv")
    df_out.to_csv(data_csv, index=False)
    print(f"[OK] Wrote data CSV: {data_csv}")

    # stats csv: x/y + label text (p and d)
    rows = []
    for i, r in enumerate(region_order):
        g1 = df[(df['Region'] == r) & (df['sim'] == 'Sim1')]['Value'].to_numpy()
        g2 = df[(df['Region'] == r) & (df['sim'] == 'Sim2')]['Value'].to_numpy()

        if len(g1) == 0 or len(g2) == 0:
            p = np.nan
            d = np.nan
            ymax = np.nan
        else:
            _, p, d = compute_wilcoxon_cohens_d(g1, g2)
            ymax = max(np.max(g1), np.max(g2))

        base = i*3 + 0.5
        x_label = base + 0.5  # centers at ticks: 1,4,7,...
        if np.isfinite(ymax):
            y_label = ymax * 1.06 + 0.01
        else:
            y_label = 0.01

        if np.isfinite(p) and np.isfinite(d):
            labtxt = f"$p={p:.2e},\\ d={d:.2f}$"
        elif np.isfinite(p) and not np.isfinite(d):
            labtxt = f"$p={p:.2e},\\ d=\\mathrm{{NA}}$"
        elif not np.isfinite(p) and np.isfinite(d):
            labtxt = f"$p=\\mathrm{{NA}},\\ d={d:.2f}$"
        else:
            labtxt = "$p=\\mathrm{NA},\\ d=\\mathrm{NA}$"

        rows.append({
            "region": r,
            "x_label": x_label,
            "y_label": y_label,
            "p_value": p,
            "cohens_d": d,
            "label_text": labtxt
        })

    stats_csv = os.path.join(output_dir, f"{metric}_stats.csv")
    pd.DataFrame(rows).to_csv(stats_csv, index=False)
    print(f"[OK] Wrote stats CSV: {stats_csv}")

# Export for every metric
for m in measure:
    export_metric_files(pd.DataFrame(all_data), m)

# Write a TeX figure per metric
for m in measure:
    write_tex_figure_for_metric(pd.DataFrame(all_data), m)

print(f"\nAll metric-specific CSVs + TeX files were created in: {output_dir}")
print("TeX files now use outlier-filtered y-range for better visualization.")