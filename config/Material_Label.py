# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:08:53 2023

@author: grife
"""

import numpy as np
import collections.abc

class Material_Label:
    """
    A class used to store data about the materials/regions to be considered

    Attributes
    ----------
    labelsMap : Map(string, arr(int))
        Dictionary mapping material name to labels used to defien that region. 
        Key: material name
        Value: array of int labels
    inverseLabelsMap : Map(int,string)
         Dictionary mnapping a single material number to the name of the material 
         Key: int material number
         Value: name of material/region

    Methods
    -------
    clear_labels_map()
        clears labelsMap and inverseLabelsMap
    removeLabel(name)
        Removes label name from labelsMap and inverseLabelsMap
    addLabelToMap(name, numbers)
        Adds material 'name' with number 'numbers' to labelsMap and inverseLabelsMap
    homogenize_material_labels(data, replace=0)
        Homogenizes data so only one label is assigned to each region. The first label in the arr givenby labelsMap
    get_homogenized_labels_map()
        Returns labels map used when homogenizing material labels in data
    
    """
    
    def __init__(self):
        self.labelsMap = {}
        self.inverseLabelsMap = {}
    
    def clear_labels_map(self):
        """
        Clears labelsMap and inverseLabelsMap

        """
        self.labelsMap = {}
        self.inverseLabelsMap = {}
    
    def removeLabel(self, name):
        """
        Removes a label name (region)

        Parameters
        ----------
        name : string
            Name of region.

        """
        name = name.lower()
        numbers = self.labelsMap.pop(name,[])
        for n in numbers:
            self.inverseLabelsMap.pop(n)

    def updateLabelInMap(self, name, numbers):

        if not (hasattr(self, "labelsMap")):
            self.create_labels_map()

        name = name.lower()

        storedLabelValues = []
        for x in list(self.labelsMap.values()):
            storedLabelValues += list(x)
        if not (isinstance(numbers, (collections.abc.Sequence, np.ndarray))):
            numbers = [numbers]

        for num in numbers:
            if storedLabelValues.count(num):
                oldName = self.inverseLabelsMap.get(num)
                self.labelsMap.get(oldName).remove(num)
                self.inverseLabelsMap.pop(num)

        self.addLabelToMap(name,numbers)


        
    def addLabelToMap(self, name, numbers):
        """
        Adds material 'name' with number 'numbers' to labelsMap and inverseLabelsMap.
        If number to added already exists, this key-value pair will not be added

        Parameters
        ----------
        name : string
            Region name.
        numbers : arr or int
            Adds elements of arr or int to labelsMap and inverseLabelsMap.

        Returns
        -------
        int
            number determining if label name was added to lists. 0 = success, otherwise failed.

        """
        if not (hasattr(self, "labelsMap")):
            self.create_labels_map()

        name = name.lower()

        storedLabelValues = []
        for x in list(self.labelsMap.values()):
            storedLabelValues += list(x)
            
        try:    
            if (isinstance(numbers, (collections.abc.Sequence, np.ndarray))):
                for num in numbers:                    
                    num = int(num)
                    if storedLabelValues.count(num):
                        print("Error. The value '{}' already exists in the map".format(num))
                        return -1
            else:
                numbers = [int(numbers)]
        except ValueError or TypeError:
            print("Error: input number '{}' not compatible".format(numbers))
            return -1
        
        labelsArr = []
        if (self.labelsMap.get(name,False)):
            labelsArr = self.labelsMap[name]
        else:
            self.labelsMap[name] = labelsArr;     
        
        labelsArr += list(numbers);
        for n in numbers:
            self.inverseLabelsMap[n] = name
        return 0;  
    
    def get_homogenized_labels_map(self):
        """
        Get labels map used when homogenizing material labels in data

        Returns
        -------
        h_label_map : Map(string, int)
            Materials mappign using only one label.

        """
        h_label_map = {}
        for key,values in self.labelsMap.items():
            h_label_map[key]=values[0]
        return h_label_map
        
    
    def homogenize_material_labels(self, data, replace=0):
        """
        Homogenizes 'data' so only one label is assigned to each region. The first label in the arr givenby labelsMap

        Parameters
        ----------
        data : array(int)
            3D array fo data to be homogenized using labelsMap property.
        replace : ine, optional
            If data valeu is not in the LabelsProperty map the homogenized value is set to this valle.
            The default is 0.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if len(data.shape) != 4:
            raise ValueError("Input data must be a 4D array.")

        print("Homogenizing 4D data according to chosen labels")
        current_dimensions = data.shape
        newData = np.zeros_like(data)
        unused_values = [replace]

        # Iterate over each voxel in the 3D space (excluding channels)
        for x in range(current_dimensions[0]):
            for y in range(current_dimensions[1]):
                for z in range(current_dimensions[2]):
                    data_value = data[x, y, z, 0]  # Access the first channel
                    
                    if self.inverseLabelsMap.get(data_value):  # Check if the value exists in the label map
                        label_name = self.inverseLabelsMap[data_value]
                        label_number = self.labelsMap[label_name][0]
                        newData[x, y, z, 0] = label_number
                    elif data_value != 0:  # Value not in label map but non-zero
                        newData[x, y, z, 0] = replace
                        if data_value not in unused_values:
                            unused_values.append(data_value)
                                    
        # Copy the second channel unchanged
        newData[..., 1] = data[..., 1]

        # Add unused values to the label map
        self.addLabelToMap("unused", unused_values)
        return newData;
            