from abc import ABC, abstractmethod
from os.path import exists
import nibabel
import os
import numpy as np

class IImport(ABC):
    @abstractmethod
    def getData(self):
        raise NotImplementedError

class ImportFromFile(IImport):
    def __init__(self, path, filename):        
        self.fullPath = "/".join([path, filename])
        assert exists(self.fullPath), f"Path {self.fullPath} does not exist."

    def getData(self):
        try:
            t1 = nibabel.load(self.fullPath)
            data = np.asarray(t1.dataobj)
            
            fa_file = os.path.join(os.getcwd(), "IOput", "in", "aseg1.mgz")
            fa_data = nibabel.load(fa_file)
            fa_values = np.asarray(fa_data.dataobj)
            print(f"data.shape: {data.shape}")
            print(f"fa_values.shape: {fa_values.shape}")
            if data.shape != fa_values.shape:
                raise ValueError("Dimensions of T1 image and FA image do not match.")

            data_with_indices = np.zeros((*data.shape, 2), dtype=np.float32)
            data_with_indices[..., 0] = data

            FA_THRESHOLD = 0.61
            high_FA_values = []
            
            def calculate_mean_of_neighbors(x, y, z):
                neighbors = []
                for dx in range(-1, 2):
                    for dy in range(-1, 2):
                        for dz in range(-1, 2):
                            if dx == 0 and dy == 0 and dz == 0:
                                continue
                            nx, ny, nz = x + dx, y + dy, z + dz
                            if 0 <= nx < fa_values.shape[0] and 0 <= ny < fa_values.shape[1] and 0 <= nz < fa_values.shape[2]:
                                neighbors.append(fa_values[nx, ny, nz] / 1000.0)
                return np.mean(neighbors) if neighbors else FA_THRESHOLD

            total_points = fa_values.size
            fa_zero_count = 0
            
            for x in range(fa_values.shape[0]):
                for y in range(fa_values.shape[1]):
                    for z in range(fa_values.shape[2]):
                        fa_value = fa_values[x, y, z] / 1000.0
                        if fa_value > FA_THRESHOLD:
                            mean_neighbor = calculate_mean_of_neighbors(x, y, z)
                            if mean_neighbor > FA_THRESHOLD:
                                high_FA_values.append((x, y, z, fa_value))
                                fa_value = FA_THRESHOLD
                            else:
                                fa_value = mean_neighbor
                        data_with_indices[x, y, z, 1] = fa_value
            
            output_file_high_fa = "./IOput/out/high_FA_values.txt"
            with open(output_file_high_fa, "w") as f:
                for value in high_FA_values:
                    f.write(f"{value}\n")
            
            output_file_data = "./IOput/out/data_with_indices.txt"
            with open(output_file_data, "w") as f:
                f.write("x y z T1 FA\n")
                for x in range(data_with_indices.shape[0]):
                    for y in range(data_with_indices.shape[1]):
                        for z in range(data_with_indices.shape[2]):
                            t1_value = data_with_indices[x, y, z, 0]
                            fa_value = data_with_indices[x, y, z, 1]
                            if t1_value != 0:
                                f.write(f"{x} {y} {z} {t1_value:.3f} {fa_value:.3f}\n")
                                if fa_value == 0:
                                    fa_zero_count += 1
            
            output_file_mgz = "./IOput/out/aseg_combined1.mgz"
            fa_values_only = data_with_indices[..., 0]
            fa_img = nibabel.Nifti1Image(fa_values_only, affine=np.eye(4))
            nibabel.save(fa_img, output_file_mgz)
            
            print("Data transfer successful.")
            print(f"Total number of points: {total_points}")
            print(f"Number of high FA values: {len(high_FA_values)}")
            print(f"Number of retained points with FA=0: {fa_zero_count}")
            
            return data_with_indices
        
        except Exception as e:
            print(f"Error while loading data: {e}")
            raise
