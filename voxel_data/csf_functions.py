import numpy as np
import save_data
from point_cloud import PointCloud
from scipy.spatial import Delaunay
from scipy import ndimage
import skimage.morphology 

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p) >= 0


class CSFFunctions:

    @staticmethod
    def add_full_csf(data, layers=1):

        current_dimensions = data.shape
        new_data = np.zeros(current_dimensions, int)

        xs, ys, zs = np.where(data == 3)
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            new_data[x, y, z] = 3

        xs, ys, zs = np.where(data == 2)
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            new_data[x, y, z] = 3

        xs, ys, zs = np.where(data == 25)
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            new_data[x, y, z] = 3

        xs, ys, zs = np.where(data == 57)
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            new_data[x, y, z] = 3

        # Create point cloud
        point_cloud = PointCloud.PointCloud()

        
        #newData_4d einfügen und mit newData vereinen.
        
        # Rufe das Array newData_4d ab
        newData_4d = save_data.get_data()
        if newData_4d is not None:
            print(f"newData_4d successfully loaded with shape: {newData_4d.shape}")
        else:
            print("No data available in newData_4d. Ensure it was saved before running this.")

        combined_array = np.zeros(new_data.shape[:3] + (2,))
        
        # copy Material value to first channel
        combined_array[..., 0] = new_data[..., -1]  # Material is lsdt dimension of new_data
        
        # Transfer FA values to the second channel only for matching coordinates
        for x in range(new_data.shape[0]):
            for y in range(new_data.shape[1]):
                for z in range(new_data.shape[2]):
                    if new_data[x, y, z] != 0:  
                        combined_array[x, y, z, 0] = new_data[x, y, z]
                        if newData_4d[x, y, z, 1] != 0:  
                            combined_array[x, y, z, 1] = newData_4d[x, y, z, 1]
                        else:  
                            combined_array[x, y, z, 1] = 0
                    else:  
                        combined_array[x, y, z, 1] = 0  
        

        new_data = combined_array



        pc = point_cloud.create_point_cloud_from_voxel(new_data)

        xmin_tot, ymin_tot, zmin_tot = [int(p) for p in np.min(pc[:, :3], axis=0)]
        xmax_tot, ymax_tot, zmax_tot = [int(p) for p in np.max(pc[:, :3], axis=0)]

        # ymax_tot = 70

        print("Filling in CSF coronal plane")
        for z in range(zmin_tot, zmax_tot + 1):
            points = point_cloud.get_slice(2, z)
            points = points[:, :2]
            if len(points) > 2:
                hull = Delaunay(points)

                xmin, ymin = np.min(points, axis=0)
                xmax, ymax_slice = np.max(points, axis=0)
                ymax = min([ymax_tot, ymax_slice])
                mid_y = int((int(ymin) + int(ymax + 1)) / 2)
                for x in range(int(xmin), int(xmax + 1)):
                    for y in range(int(ymin), int(ymax + 1)):
                        if data[x, y, z] == 0 and (y < ymax_tot):
                            if in_hull([x, y], hull):
                                # add FA-Value (0)
                                point_cloud.add_point_to_cloud([x, y, z, 24, 0])
                                data[x, y, z] = 24
                                new_data[x, y, z] = 24
            else:
                print("y-value", z)
                print("points: ", points)

        ##### CSF is not filled for the sagital plane of the brain
        # print("Filling in CSF x-dim")
        # for x in range(xmin_tot, xmax_tot + 1):
        #     points = point_cloud.get_slice(0, x)
        #     points = points[:, 1:3]
        #     hull = Delaunay(points)
        #
        #     min1d, min2d = np.min(points, axis=0)
        #     max1d_slice, max2d = np.max(points, axis=0)
        #     max1d = min([ymax_tot, max1d_slice])
        #     mid_2d = int((int(min2d) + int(max2d + 1)) / 2)
        #     for y in range(int(min1d), int(max1d + 1)):
        #         for z in range(int(min2d), int(max2d + 1)):
        #             if (data[x, y, z] == 0) and (y < ymax_tot):
        #                 if in_hull([y, z], hull):
        #                     point_cloud.add_point_to_cloud([x, y, z, 24])
        #                     data[x, y, z] = 24
        #                     new_data[x, y, z] = 24

        print("Filling in CSF transverse plane")
        for y in range(ymin_tot, ymax_tot + 1):
            points = point_cloud.get_slice(1, y)
            points = points[:, [0, 2]]
            if len(points) > 2:
                hull = Delaunay(points)

                min1d, min2d = np.min(points, axis=0)
                max1d, max2d = np.max(points, axis=0)
                mid_2d = int((int(min2d) + int(max2d + 1)) / 2)
                for x in range(int(min1d), int(max1d + 1)):
                    for z in range(int(min2d), int(max2d + 1)):
                        if (data[x, y, z] == 0) and (y < ymax_tot):
                            if in_hull([x, z], hull):
                                # add FA-Value (0)
                                point_cloud.add_point_to_cloud([x, y, z, 24, 0])
                                data[x, y, z] = 24
                                new_data[x, y, z] = 24
            else:
                print("y-value", y)
                print("points: ", points)


##### Old method - used by benedikt. compare to results of new method below.
        # inflated_csf = ndimage.binary_dilation(data).astype(int)
        # for i in range(layers - 1):
        #     inflated_csf = ndimage.binary_dilation(inflated_csf).astype(int)

        # xs, ys, zs = np.where(inflated_csf == 1)
        # current_dimensions = inflated_csf.shape
        # for [x, y, z] in np.column_stack((xs, ys, zs)):
        #     if data[x, y, z] == 0:
        #         data[x, y, z] = 24


        structure_mask = (data > 0).astype(bool)
        
        # Method 1: Heavy Gaussian pre-smoothing of the mask
        smoothed_mask = ndimage.gaussian_filter(structure_mask.astype(float), sigma=3.0)
        distance_from_structure = ndimage.distance_transform_edt(smoothed_mask < 0.5)
        inflated_csf = distance_from_structure <= layers
        
        # Method 2: Multiple rounds of heavy smoothing
        for _ in range(3):
            inflated_csf = ndimage.gaussian_filter(inflated_csf.astype(float), sigma=2.0) > 0.3
        
        # Method 3: Final morphological smoothing
        sphere = ndimage.generate_binary_structure(3, 2)  # Higher connectivity
        inflated_csf = ndimage.binary_closing(inflated_csf, structure=sphere, iterations=3)
        inflated_csf = ndimage.binary_opening(inflated_csf, structure=sphere, iterations=2)
        
        # Alternative ultra-smooth method (uncomment to try instead):
        # # Convolution with large smooth kernel
        # kernel_size = 7
        # kernel = np.ones((kernel_size, kernel_size, kernel_size)) / (kernel_size**3)
        # smoothed = ndimage.convolve(structure_mask.astype(float), kernel)
        # distance_from_structure = ndimage.distance_transform_edt(smoothed < 0.3)
        # inflated_csf = distance_from_structure <= layers
        # inflated_csf = ndimage.gaussian_filter(inflated_csf.astype(float), sigma=3.0) > 0.2
        
        xs, ys, zs = np.where(inflated_csf == 1)
        current_dimensions = inflated_csf.shape
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            if data[x, y, z] == 0:
                data[x, y, z] = 24

    @staticmethod
    def add_partial_csf(data, layers=1):

        current_dimensions = data.shape
        new_data = np.zeros(current_dimensions, int)

        xs, ys, zs = np.where(data == 3)
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            new_data[x, y, z] = 24

        xs, ys, zs = np.where(data == 2)
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            new_data[x, y, z] = 24

        xs, ys, zs = np.where(data == 25)
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            new_data[x, y, z] = 24

        xs, ys, zs = np.where(data == 57)
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            new_data[x, y, z] = 24

        # xs, ys, zs = np.where(new_data == 1)
        # current_dimensions = new_data.shape
        # for [x, y, z] in np.column_stack((xs, ys, zs)):
        #     if new_data[x, y, z] == 0:
        #         new_data[x, y, z] = 24

                # Create point cloud
        point_cloud = PointCloud.PointCloud()
        pc = point_cloud.create_point_cloud_from_voxel(new_data)

        xmin_tot, ymin_tot, zmin_tot = [int(p) for p in np.min(pc[:, :3], axis=0)]
        xmax_tot, ymax_tot, zmax_tot = [int(p) for p in np.max(pc[:, :3], axis=0)]

        print("Filling in CSF coronal plane")
        for z in range(zmin_tot, zmax_tot + 1):
            points = point_cloud.get_slice(2, z)
            points = points[:, :2]
            if len(points) > 2:
                hull = Delaunay(points)

                xmin, ymin = np.min(points, axis=0)
                xmax, ymax_slice = np.max(points, axis=0)
                ymax = min([ymax_tot, ymax_slice])
                mid_y = int((int(ymin) + int(ymax + 1)) / 2)
                for x in range(int(xmin), int(xmax + 1)):
                    for y in range(int(ymin), int(ymax + 1)):
                        if (data[x, y, z] == 0) and (y < ymax_tot):
                            if in_hull([x, y], hull):
                                new_data[x, y, z] = 24
            else:
                print("z-value", z)
                print("points: ", points)


        ##### CSF is not filled for the sagital plane of the brain
        # print("Filling in CSF x-dim")
        # for x in range(xmin_tot, xmax_tot + 1):
        #     points = point_cloud.get_slice(0, x)
        #     points = points[:, 1:3]
        #     hull = Delaunay(points)
        #
        #     min1d, min2d = np.min(points, axis=0)
        #     max1d_slice, max2d = np.max(points, axis=0)
        #     max1d = min([ymax_tot, max1d_slice])
        #     mid_2d = int((int(min2d) + int(max2d + 1)) / 2)
        #     for y in range(int(min1d), int(max1d + 1)):
        #         for z in range(int(min2d), int(max2d + 1)):
        #             if (data[x, y, z] == 0) and (y < ymax_tot):
        #                 if in_hull([y, z], hull):
        #                     new_data[x, y, z] = 24

        print("Filling in CSF trasverse plane")
        for y in range(ymin_tot, ymax_tot + 1):
            points = point_cloud.get_slice(1, y)
            points = points[:, [0, 2]]
            if len(points) > 2:
                hull = Delaunay(points)

                min1d, min2d = np.min(points, axis=0)
                max1d, max2d = np.max(points, axis=0)
                mid_2d = int((int(min2d) + int(max2d + 1)) / 2)
                for x in range(int(min1d), int(max1d + 1)):
                    for z in range(int(min2d), int(max2d + 1)):
                        if (data[x, y, z] == 0) and (y < ymax_tot):
                            if in_hull([x, z], hull):
                                new_data[x, y, z] = 24
            else:
                print("y-value", y)
                print("points: ", points)

        new_data = ndimage.binary_erosion(new_data).astype(int)
        # for i in range(1):
        #     new_data = ndimage.binary_erosion(new_data).astype(int)

        xs, ys, zs = np.where(new_data == 1)
        current_dimensions = new_data.shape
        for [x, y, z] in np.column_stack((xs, ys, zs)):
            if data[x, y, z] == 0:
                data[x, y, z] = 24

    @staticmethod
    def get_csf_function(csf_type):
        if csf_type == "full":
            return CSFFunctions.add_full_csf
        elif csf_type == "partial":
            return CSFFunctions.add_partial_csf
        else:
            return None
