import os
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
matplotlib.use("Agg")       
import matplotlib.pyplot as plt
import pandas as pd
import csv



INPUT_DIRECTORY_PATH = './output_atrophy/post_processing/rampp_final/brainrampp/'
OUTPUT_DIRECTORY_PATH = './output_atrophy/post_processing/rampp_final/plots/volume_evolution_per_region/'

BASE_FILENAMES = ['rampp9R', 'ramppFA']  # models

# Brain material IDs
BRAIN_MATERIAL_IDS = [251, 1024, 13, 2, 12, 1011, 7, 173, 10, 11, 18, 174, 17, 175, 1035, 1028, 1030]

# Ensure output directory exists
os.makedirs(OUTPUT_DIRECTORY_PATH, exist_ok=True)

# Function to compute volume of a hexahedral cell
def compute_hex_volume(coords):
    if len(coords) == 8:  # Ensure hexahedral cell
        v0 = coords[0]
        return abs(np.dot(np.cross(coords[1] - v0, coords[2] - v0), coords[4] - v0)) / 6.0
    return 0.0


LABELS = {
    'rampp9R': '9R model',
    'ramppFA': 'FA model',
}

# Initialize lists to store results for all base filenames
all_steps = []
all_brain_volumes = []
all_brain_volume_changes = []
# all_volume_fractions = []

# Loop over each BASE_FILENAME
for BASE_FILENAME in BASE_FILENAMES:
    # Initialize storage for this base filename
    steps = []
    brain_volumes = []
    brain_volume_changes = []
    # volume_fractions = []


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
        material_volumes = {material_id: 0.0 for material_id in BRAIN_MATERIAL_IDS}
        material_volumes[24] = 0.0  # CSF
        material_volumes[4] = 0.0  # Ventricles
        material_volumes['entire_volume'] = 0.0

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
            if material_id in BRAIN_MATERIAL_IDS:
                material_volumes[material_id] += cell_volume
            elif material_id == 24:  # Exclude CSF
                material_volumes[24] += cell_volume
            elif material_id == 4:  # Exclude ventricles
                material_volumes[4] += cell_volume

        # Adjust brain volume by excluding material ID 24 (CSF) and 4 (ventricles)
        brain_volume = material_volumes['entire_volume'] - material_volumes[24] - material_volumes[4]
        print(f"Brain volume (step {step}): {brain_volume}")

        # Store computed values
        brain_volumes.append(brain_volume)
        if step == 0:
            initial_brain_volume = brain_volume
        brain_volume_changes.append(100 * (brain_volume - initial_brain_volume) / initial_brain_volume)
        # volume_fractions.append(100 * brain_volume / material_volumes['entire_volume'])
        steps.append(step * 5)  # Modify x-axis to be 0, 5, 10, ..., 80

    # Append results for this base filename to the overall lists
    all_steps.append(steps)
    all_brain_volumes.append(brain_volumes)
    all_brain_volume_changes.append(brain_volume_changes)
    # all_volume_fractions.append(volume_fractions)

# Create combined plots
for i, (data, ylabel, title, ylims) in enumerate([
    (all_brain_volumes, "Brain Volume", "Brain Volume Evolution", None),
    (all_brain_volume_changes, "Percentage Change (%)", "Brain Volume Percentage Change", None),
    # (all_volume_fractions, "Volume Fraction Percentage (%)", "Brain Volume Fraction Percentage", None)
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

# Export brain volume percentage change data to tables
print("\n=== Exporting Brain Volume Percentage Change Tables ===")

# Method 1: Export to CSV
csv_output_path = os.path.join(OUTPUT_DIRECTORY_PATH, "brain_volume_percentage_change.csv")
with open(csv_output_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    
    # Write header
    header = ['Step'] + [LABELS[base_filename] for base_filename in BASE_FILENAMES]
    writer.writerow(header)
    
    # Write data rows
    for i in range(len(all_steps[0])):  # Assuming all models have the same number of steps
        row = [all_steps[0][i]]  # Step value
        for j in range(len(BASE_FILENAMES)):
            row.append(f"{all_brain_volume_changes[j][i]:.4f}")
        writer.writerow(row)

print(f"Brain volume percentage change data exported to CSV: {csv_output_path}")

# Method 2: Export to pandas DataFrame and save as Excel (if pandas is available)
try:
    # Create DataFrame
    df_dict = {'Step': all_steps[0]}
    for i, base_filename in enumerate(BASE_FILENAMES):
        df_dict[LABELS[base_filename]] = all_brain_volume_changes[i]
    
    df = pd.DataFrame(df_dict)
    
    # Save to Excel
    excel_output_path = os.path.join(OUTPUT_DIRECTORY_PATH, "brain_volume_percentage_change.xlsx")
    df.to_excel(excel_output_path, index=False, sheet_name='Brain Volume Change')
    print(f"Brain volume percentage change data exported to Excel: {excel_output_path}")
    
    # Display table summary
    print("\n=== Brain Volume Percentage Change Summary ===")
    print(df.round(4))
    
except ImportError:
    print("Pandas not available. Excel export skipped. CSV export completed successfully.")

# Method 3: Export formatted text table
txt_output_path = os.path.join(OUTPUT_DIRECTORY_PATH, "brain_volume_percentage_change.txt")
with open(txt_output_path, 'w') as txtfile:
    txtfile.write("Brain Volume Percentage Change Table\n")
    txtfile.write("="*50 + "\n\n")
    
    # Write header
    header = f"{'Step':<6}"
    for base_filename in BASE_FILENAMES:
        header += f"{LABELS[base_filename]:<15}"
    txtfile.write(header + "\n")
    txtfile.write("-" * len(header) + "\n")
    
    # Write data rows
    for i in range(len(all_steps[0])):
        row = f"{all_steps[0][i]:<6}"
        for j in range(len(BASE_FILENAMES)):
            row += f"{all_brain_volume_changes[j][i]:<15.4f}"
        txtfile.write(row + "\n")
    
    # Add summary statistics
    txtfile.write("\n" + "="*50 + "\n")
    txtfile.write("SUMMARY STATISTICS\n")
    txtfile.write("="*50 + "\n\n")
    
    for i, base_filename in enumerate(BASE_FILENAMES):
        changes = all_brain_volume_changes[i]
        txtfile.write(f"{LABELS[base_filename]}:\n")
        txtfile.write(f"  Initial Change: {changes[0]:.4f}%\n")
        txtfile.write(f"  Final Change:   {changes[-1]:.4f}%\n")
        txtfile.write(f"  Max Change:     {max(changes):.4f}%\n")
        txtfile.write(f"  Min Change:     {min(changes):.4f}%\n")
        txtfile.write(f"  Average Change: {np.mean(changes):.4f}%\n\n")

print(f"Brain volume percentage change data exported to formatted text: {txt_output_path}")

print("\n=== Export Summary ===")
print(f"All brain volume percentage change data has been exported to:")
print(f"1. CSV format: {csv_output_path}")
print(f"2. Formatted text: {txt_output_path}")
if 'excel_output_path' in locals():
    print(f"3. Excel format: {excel_output_path}")
print("\nTable export completed successfully!")