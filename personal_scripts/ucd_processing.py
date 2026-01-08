import os
import pandas as pd

def process_UCD_extension(input_file, output_file):
    # Read the entire input file
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Extract header lines (starting with "#")
    header_lines = [line for line in lines if line.startswith("#")]
    # Get the remaining lines (without header)
    non_header_lines = [line for line in lines if not line.startswith("#")]

    # Check if the first non-header line is a metadata line:
    metadata_line = None
    if non_header_lines:
        tokens = non_header_lines[0].split()
        if len(tokens) == 5 and tokens[0] != "0":
            metadata_line = non_header_lines.pop(0)

    # Determine from which line the element data starts (i.e., the 3rd token is "hex" or "quad")
    element_start = next(i for i, line in enumerate(non_header_lines)
                          if len(line.split()) >= 3 and line.split()[2] in ["hex", "quad"])
    node_lines = non_header_lines[:element_start]
    element_lines = non_header_lines[element_start:]

    # Process node data: If a line has only 4 tokens, append an FA value "0"
    processed_node_lines = []
    for line in node_lines:
        tokens = line.split()
        if len(tokens) == 4:
            tokens.append("0")
        if len(tokens) == 5:
            processed_node_lines.append("\t".join(tokens) + "\n")
        else:
            continue

    # Load node data into a DataFrame
    node_data = [line.strip().split() for line in processed_node_lines]
    node_df = pd.DataFrame(node_data, columns=["Nummer", "x", "y", "z", "FA"]).astype({
        "Nummer": int, "x": float, "y": float, "z": float, "FA": float
    })

    def calculate_mean_fa(element_nodes):
        # Consider only nodes with FA not equal to 0
        fa_values = node_df[node_df["Nummer"].isin(element_nodes) & (node_df["FA"] != 0)]["FA"]
        return fa_values.mean() if not fa_values.empty else 0

    # Process element data for hex and quad elements
    element_list = []
    for line in element_lines:
        tokens = line.split()
        if len(tokens) < 3:
            continue
        element_type = tokens[2]
        if element_type == "hex" and len(tokens) >= 11 and all(tokens[i].isdigit() for i in [0, 1] + list(range(3, 11))):
            nodes = [int(token) for token in tokens[3:11]]
            mean_fa = calculate_mean_fa(nodes)
            element_list.append((tokens, mean_fa))
        elif element_type == "quad" and len(tokens) >= 7 and all(tokens[i].isdigit() for i in [0, 1] + list(range(3, 7))):
            nodes = [int(token) for token in tokens[3:7]]
            mean_fa = calculate_mean_fa(nodes)
            element_list.append((tokens, mean_fa))

    # Write the updated file
    with open(output_file, "w") as f:
        for line in header_lines:
            f.write(line)
        if metadata_line:
            f.write(metadata_line)
        for line in processed_node_lines:
            f.write(line)
        for tokens, mean_fa in element_list:
            element_line = "\t".join(tokens) + f"\t{mean_fa:.6f}\n"
            f.write(element_line)

    print(f"Putput file: '{output_file}'.")


def process_UCD_column_remove(input_file, output_file):

#    Removes the last column in lines with exactly 5 columns (split by tabs),   except for the first data line.

    first_data_line = True  # The first data line (without header/comments) is copied unchanged
    with open(input_file, 'r', encoding='utf-8') as f_in, open(output_file, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            # Copy empty lines directly
            if not line.strip():
                f_out.write(line)
                continue
            # Copy comment lines
            if line.strip().startswith('#'):
                f_out.write(line)
                continue
            # Copy the first data line completely
            if first_data_line:
                f_out.write(line)
                first_data_line = False
                continue
            # For lines with exactly 5 columns (split by tabs), remove the last column
            tokens = line.rstrip('\n').split('\t')
            if len(tokens) == 5:
                new_line = '\t'.join(tokens[:-1]) + '\n'
                f_out.write(new_line)
            else:
                f_out.write(line)
    print(f"Output file:'{output_file}'.")


def convert_to_vtk(input_filename, output_filename):
    nodes = {}
    elements = []
    cell_scalars = []

    with open(input_filename, 'r') as f:
        for line in f:
            tokens = line.strip().split()
            if not tokens:
                continue

            if len(tokens) == 4:
                try:
                    node_id = int(tokens[0])
                    x, y, z = map(float, tokens[1:4])
                    nodes[node_id] = [x, y, z]
                except ValueError:
                    continue

            elif len(tokens) == 12 and tokens[2].lower() == "hex":
                try:
                    elem_id = int(tokens[0])
                    region_id = int(tokens[1])
                    node_ids = list(map(int, tokens[3:11]))
                    scalar_val = float(tokens[11])
                    elements.append({
                        "elem_id": elem_id,
                        "region_id": region_id,
                        "type": "hex",
                        "node_ids": node_ids,
                        "scalar": scalar_val
                    })
                    cell_scalars.append(scalar_val)
                except ValueError:
                    continue

            elif len(tokens) == 8 and tokens[2].lower() == "quad":
                try:
                    elem_id = int(tokens[0])
                    region_id = int(tokens[1])
                    node_ids = list(map(int, tokens[3:7]))
                    scalar_val = float(tokens[7])
                    elements.append({
                        "elem_id": elem_id,
                        "region_id": region_id,
                        "type": "quad",
                        "node_ids": node_ids,
                        "scalar": scalar_val
                    })
                    cell_scalars.append(scalar_val)
                except ValueError:
                    continue
            else:
                continue

    if not nodes:
        print("Keine Knotendaten gefunden!")
        exit(1)
    if not elements:
        print("Keine Elementdaten gefunden!")
        exit(1)

    sorted_node_ids = sorted(nodes.keys())
    node_to_index = {nid: idx for idx, nid in enumerate(sorted_node_ids)}

    with open(output_filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Converted mesh with cell scalars\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")

        num_points = len(sorted_node_ids)
        f.write(f"POINTS {num_points} float\n")
        for nid in sorted_node_ids:
            x, y, z = nodes[nid]
            f.write(f"{x} {y} {z}\n")

        num_cells = len(elements)
        total_ints = sum([(len(elem["node_ids"]) + 1) for elem in elements])
        f.write(f"\nCELLS {num_cells} {total_ints}\n")
        for elem in elements:
            n_nodes = len(elem["node_ids"])
            indices = [node_to_index[nid] for nid in elem["node_ids"]]
            f.write(f"{n_nodes} " + " ".join(str(idx) for idx in indices) + "\n")

        f.write(f"\nCELL_TYPES {num_cells}\n")
        for elem in elements:
            if elem["type"].lower() == "hex":
                f.write("12\n")  # VTK_HEX_Element
            elif elem["type"].lower() == "quad":
                f.write("9\n")   # VTK_QUAD_Element

        f.write(f"\nCELL_DATA {num_cells}\n")
        f.write("SCALARS cell_scalar float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for scalar in cell_scalars:
            f.write(f"{scalar}\n")

    print("output file:", output_filename)

