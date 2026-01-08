from readers.Readers import Reader
from writers.Writers import Writer
from mesh.Node import Node
from mesh.Mesh import Mesh
from mesh.Element import Element
from vtk import vtkIdList, vtkUnstructuredGridReader
from os.path import exists
from mesh.mesh_transformations import translate_mesh


### This scripts transforms the data to show the deformation in the current configuration. 
### Takes files from extract_last_step or manually extracted files

def readVtk(path, filename):
    fullFilename1 = path + "/" + filename + ".vtk"
    if (exists(fullFilename1)):
        print("##Reading: " + filename)

        # Set up poly data reader for result set 1
        reader = vtkUnstructuredGridReader()
        reader.SetFileName(fullFilename1)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        grid = reader.GetOutput()
        return grid
    print("Could not find file with name {}".format(fullFilename1))
    raise FileNotFoundError

def create_deformation_plot(f, path):
    inputPath = path
    dataMap = {}
    grid = readVtk(inputPath, f)

    mesh = Mesh()
    mesh.dataToWrite = ["displacement"]
    mesh.cellData = ["displacement_centroid"]
    displacementArray = grid.GetPointData().GetArray("displacement")
    materialsArray = grid.GetPointData().GetArray("material_ids")

    disp_data = {}  # Key = starting position, displacement
    point_position_to_nodes = {}  # Key= point position, value  = node value
    cell_point_to_nodes = {}  # Key= cell pointId, value  = node value
    nodeMap = {}
    elementsMap = {}
    elementToMaterial = {}
    for k in range(grid.GetNumberOfCells()):
        cellIds = vtkIdList()
        grid.GetCellPoints(k, cellIds)
        ica = []
        mat = 0
        displacement_tot = [0, 0, 0];
        # von_Mises = 0
        for i in range(8):
            pointId = cellIds.GetId(i)
            pointPosition = grid.GetPoint(pointId)
            mat += materialsArray.GetValue(pointId)
            position_key = "-".join([str(x) for x in pointPosition])
            nodeValue = point_position_to_nodes.get(position_key, -1)
            if nodeValue != -1:
                cell_point_to_nodes[pointId] = nodeValue
                node = nodeMap[int(nodeValue)]
                displacement = disp_data[nodeValue];
                assert [disp_data.get(nodeValue, -1) != -1]
            else:
                nodeValue = pointId
                cell_point_to_nodes[pointId] = nodeValue
                point_position_to_nodes[position_key] = pointId
                displacement = [round(y, 6) for y in displacementArray.GetTuple(pointId)]
                # von_Mises += round(vonMisesArray.GetTuple(pointId)[0], 6)
                disp_data[pointId] = displacement
                newPosition = [round(y, 6) for y in pointPosition]
                for d in range(len(displacement)):
                    newPosition[d] = pointPosition[d] + displacement[d]
                    if d == 2:
                        newPosition[d] += 4
                node = Node(int(nodeValue), newPosition)
                node.addData("displacement", displacement)
                nodeMap[int(nodeValue)] = node
            for d in range(3):
                displacement_tot[d] += displacement[d]
            ica.append(node)
        element = Element(k, ica)
        element.setMaterial((mat / 8))
        element.properties['displacement_centroid'] = [d / 8. for d in displacement_tot]
        # element.properties['von_mises'] = [von_Mises / 8.]
        elementsMap[k] = element

    mesh.nodes = nodeMap
    mesh.elements = elementsMap

    translate_mesh(mesh.nodes, [0,0,-4])

    ## Write deformed mesh to new vtk
    path = inputPath
    filename = f + "_deformed"
    writer = Writer()
    writer.openWriter("vtk", filename, path)
    writer.writeMeshData(mesh)
    writer.closeWriter()

if __name__ == "__main__":    
    path = "./output_atrophy/post_processing/rampp_FA_9R/"
    # path = "./output_atrophy/OAS1_0002_MR1_atrophy_fine"
    # filenames = ["OAS1_0002_MR1_atrophy_1R_last_timestep", "OAS1_0002_MR1_atrophy_2R_last_timestep", "OAS1_0002_MR1_atrophy_4R_last_timestep", "OAS1_0002_MR1_atrophy_9R_last_timestep", "OAS1_0002_MR1_atrophy_19R_last_timestep"]
    filenames = ["rampp_9R_last_timestep","rampp_FA_last_timestep"]
    for f in filenames:
        create_deformation_plot(f, path)
