class ArrayProcessor:
    def __init__(self):
        pass

    def printArray(self, array, word):

        import numpy as np
        import os


        script_dir = os.path.dirname(os.path.abspath(__file__))
        relative_path = os.path.join("IOput", "out", "atrophy_files", "rampp")
        output_path = os.path.join(script_dir, relative_path)

        # Filepath
        #output_path = r"C:\Users\pumab\Documents\Studium Master\Masterarbeit\3_Implementation_dealII\Brain_Creation_Code-17-region-model\TestPrint"  # Das Verzeichnis, in dem die Datei gespeichert wird

        if not isinstance(array, np.ndarray):
            raise ValueError("Array must be a numpy-Array.")

        if array.ndim not in [3, 4]:
            raise ValueError("Array has to be 3D or 4D.")

        output_lines = []
        shape = array.shape

        non_zero_count = 0  # Zähler für Punkte ungleich 0
        total_count = np.prod(shape[:3])  

        if array.ndim == 3:
            for z in range(shape[0]):
                for y in range(shape[1]):
                    for x_coord in range(shape[2]):
                        value = array[z, y, x_coord]
                        if value != 0:
                            non_zero_count += 1
                            output_lines.append(f"{x_coord} {y} {z} {value}")
        elif array.ndim == 4:
            for z in range(shape[0]):
                for y in range(shape[1]):
                    for x_coord in range(shape[2]):
                        material_value = array[z, y, x_coord, 0]
                        fa_value = array[z, y, x_coord, 1]
                        if material_value != 0:
                            non_zero_count += 1
                            output_lines.append(
                                f"{x_coord} {y} {z} {material_value} {fa_value}"
                            )

        if not os.path.exists(output_path):
            os.makedirs(output_path)

        output_file = os.path.join(output_path, f"{word}_Array.txt")

        with open(output_file, "w") as file:
            file.write(f"Anzahl der Punkte ungleich 0: {non_zero_count}\n")
            file.write(f"Gesamtanzahl der Punkte: {total_count}\n\n")
            for line in output_lines:
                file.write(line + "\n")

        print(f"Das Array wurde in der Datei {output_file} gespeichert.")

