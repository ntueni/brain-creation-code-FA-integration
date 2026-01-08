import os
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

# Force headless backend BEFORE importing pyplot
import matplotlib
matplotlib.use("Agg")   # or set env var MPLBACKEND=Agg before running
import matplotlib.pyplot as plt


# Input and output paths
INPUT_DIRECTORY_PATH = './output_atrophy/post_processing/rampp_FA_9R/brainrampp/'
OUTPUT_DIRECTORY_PATH = './output_atrophy/post_processing/rampp_FA_9R/plots/'


BASE_FILENAMES = [
    'rampp9R',
    'ramppFA'
]

TARGET_MATERIAL_ID = 24  # CSF

# Ensure output directory exists
os.makedirs(OUTPUT_DIRECTORY_PATH, exist_ok=True)

# Function to compute volume of a hexahedral cell
def compute_hex_volume(coords):
    if len(coords) == 8:  # Ensure hexahedral cell
        v0 = coords[0]
        return abs(np.dot(np.cross(coords[1] - v0, coords[2] - v0), coords[4] - v0)) / 6.0
    return 0.0

# Mapping base filenames to readable labels
# LABELS = {
#     'OAS1_0002_MR1_atrophy_1R': '1R model',
#     'OAS1_0002_MR1_atrophy_2R': '2R model',
#     'OAS1_0002_MR1_atrophy_4R': '4R model',
#     'OAS1_0002_MR1_atrophy_9R': '9R model',
#     'OAS1_0002_MR1_atrophy_19R': '19R model'
# }

LABELS = {
    'rampp9R' : '9R model',
    'rampp17' : '17R model',
    'ramppFA' : 'FA model'
}

# Initialize lists to store results for all base filenames
all_steps = []
all_CSF_volumes = []
all_CSF_volume_changes = []
all_volume_fractions = []

# Loop over each BASE_FILENAME
for BASE_FILENAME in BASE_FILENAMES:
    # Initialize storage for this base filename
    steps = []
    CSF_volumes = []
    CSF_volume_changes = []
    volume_fractions = []
    brain_volumes = []

    # Loop through the steps
    for step in range(17):
        vtk_file = os.path.join(INPUT_DIRECTORY_PATH, f"{BASE_FILENAME}_{step}.vtk")
        print(f"Processing {vtk_file}...")

        # Read the VTK file
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtk_file)
        reader.Update()
        data = reader.GetOutput()

        # Find material IDs array
        material_ids_array = None
        for idx in range(data.GetPointData().GetNumberOfArrays()):
            array = data.GetPointData().GetArray(idx)
            if "material" in array.GetName().lower():
                material_ids_array = array
                break

        if material_ids_array is None:
            raise ValueError(f"No material IDs array found in {vtk_file}")

        # Convert material IDs to NumPy array
        material_ids_np = vtk_to_numpy(material_ids_array)

        # Initialize volumes
        material_volumes = {TARGET_MATERIAL_ID: 0.0, 24: 0.0, 4: 0.0, 'entire_volume': 0.0}

        # Process each cell
        for cell_id in range(data.GetNumberOfCells()):
            cell = data.GetCell(cell_id)

            # Get the material ID for the cell (first point's ID as representative)
            material_id = int(material_ids_np[cell.GetPointIds().GetId(0)])

            # Get cell points and compute volume
            points = cell.GetPoints()
            coords = np.array([points.GetPoint(i) for i in range(points.GetNumberOfPoints())])

            cell_volume = compute_hex_volume(coords)

            # Update total volume and material-specific volumes
            material_volumes['entire_volume'] += cell_volume
            if material_id == TARGET_MATERIAL_ID:
                material_volumes[TARGET_MATERIAL_ID] += cell_volume
            # elif material_id in [24, 4]:  # Exclude material IDs 24 (CSF) and 4 (ventricles)
            #     material_volumes[material_id] += cell_volume

        # Adjust brain volume by excluding material IDs 24 and 4
        brain_volume = material_volumes['entire_volume'] - material_volumes[24] - material_volumes[4]
        print("brain_volume = ")
        print(brain_volume)

        # Store computed values
        CSF_volume = material_volumes[TARGET_MATERIAL_ID]
        CSF_volumes.append(CSF_volume)
        volume_fractions.append(100*CSF_volume / brain_volume)
        if step == 0:
            initial_CSF_volume = CSF_volume
        CSF_volume_changes.append(100 * (CSF_volume-initial_CSF_volume) / initial_CSF_volume)
        steps.append(step * 5)  # Modify x-axis to be 0, 5, 10, ..., 80

    # Append results for this base filename to the overall lists
    all_steps.append(steps)
    all_CSF_volumes.append(CSF_volumes)
    all_CSF_volume_changes.append(CSF_volume_changes)
    all_volume_fractions.append(volume_fractions)

# Create combined plots
for i, (data, ylabel, title, ylims) in enumerate([
    (all_CSF_volumes, "CSF Volume", "CSF Volume Evolution", None),
    (all_CSF_volume_changes, "CSF Change (%)", "CSF Volume Percentage Change", None),
    (all_volume_fractions, "Volume Fraction Percentage (%)", "CSF Volume Fraction Percentage", None)
]):
    plt.figure(figsize=(12, 6))

    # Loop through each base filename's data and plot it
    for j, base_filename in enumerate(BASE_FILENAMES):
        plt.plot(all_steps[j], data[j], marker='o', label=LABELS[base_filename])

    # Customize the plot
    plt.xlabel("Steps")
    plt.ylabel(ylabel)
    plt.title(f"{title}")
    plt.grid(True, which='both', axis='both', linestyle='--', color='gray') 
    plt.legend(title="Model")

    # Set the y-axis limits for percentage change plots
    if ylims:
        plt.ylim(ylims)

    # Save the plot as PNG
    output_plot_path = os.path.join(OUTPUT_DIRECTORY_PATH, f"{title.replace(' ', '_')}.png")
    plt.tight_layout()
    plt.savefig(output_plot_path)
    plt.close()  # Close the plot to save memory

    print(f"Plot saved as {output_plot_path}")