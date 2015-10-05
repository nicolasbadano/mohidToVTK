#include "MohidResults.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include "InputFile.h"

MohidResults::MohidResults(double offset)
{
    verticalOffset = offset;
    mapIsLoaded = false;

    fields.reserve(5);

    hasSalinity = false;
    hasDensity = false;
    hasColiforms = false;
    hasTemperature = false;
}

MohidResults::~MohidResults(void)
{
}

//---------------------------------------------------------------------------------
//----------- loadHydrodynamicResults ---------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadHydrodynamicResults(char* hdf5FileName, int index)
{
    int col, row, layer;

    hid_t       file_id;        // File handle
    herr_t      status;         // Status of the functions
    hsize_t     dims[3];        // Dimensions of the array
    float       *data;
    int         *maskData;
    char        datasetName[200];

    cout << "\n";
    getDatasetName("", index, datasetName);
    cout << "Step " << datasetName << "- File " << hdf5FileName << " - Loading hydrodynamic results..." << "\n";

    // Open HDF5 file
    file_id = H5Fopen (hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Get mask dimensions
    getDatasetName("/Grid/OpenPoints/OpenPoints_", index, datasetName);
    status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
    if (status == -1) return false;
    numLayers = dims[0];
    numCol = dims[1];
    numRow = dims[2];

    maskData = (int *)malloc(sizeof(int)*(numCol+1)*(numRow+1)*(numLayers+1));
    data = (float *)malloc(sizeof(float)*(numCol+1)*(numRow+1)*(numLayers+1));

    // Read mask data
    status = H5LTread_dataset_int(file_id,datasetName, maskData);

    // Store the mask data
    mask.resize(numLayers);
    for (layer = 0; layer < numLayers; layer++) {
        mask[layer].resize(numCol);
        for (col = 0; col < numCol; col++) {
            mask[layer][col].resize(numRow,(int)0);
            for (row = 0; row < numRow; row++) {
                mask[layer][col][row] = maskData[layer*numCol*numRow+col*numRow+row];
            }
        }
    }

    // Get ConnectionX dimensions
    status = H5LTget_dataset_info(file_id,"/Grid/ConnectionX",dims,NULL,NULL);
    if (status == -1) return false;
    if (dims[0] != numCol+1) return false;
    if (dims[1] != numRow+1) return false;

    // Read ConnectionX dataset
    status = H5LTread_dataset_float(file_id,"/Grid/ConnectionX", data);

    x.resize(numCol+1);
    for (col = 0; col < numCol+1; col++) {
        x[col].resize(numRow+1,(float)0);
        for (row = 0; row < numRow+1; row++) {
            x[col][row] = data[col*(numRow+1)+row];
        }
    }

    // Get ConnectionY dimensions
    status = H5LTget_dataset_info(file_id,"/Grid/ConnectionY",dims,NULL,NULL);
    if (status == -1) return false;
    if (dims[0] != numCol+1) return false;
    if (dims[1] != numRow+1) return false;

    // Read ConnectionY dataset
    status = H5LTread_dataset_float(file_id,"/Grid/ConnectionY", data);

    y.resize(numCol+1);
    for (col = 0; col < numCol+1; col++) {
        y[col].resize(numRow+1,(float)0);
        for (row = 0; row < numRow+1; row++) {
            y[col][row] = data[col*(numRow+1)+row];
        }
    }


    // Get VerticalZ dimensions
    getDatasetName("/Grid/VerticalZ/Vertical_", index, datasetName);
    status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
    if (status == -1) return false;
    if (dims[0] != numLayers+1) return false;
    if (dims[1] != numCol) return false;
    if (dims[2] != numRow) return false;

    // Read VerticalZ dataset
    status = H5LTread_dataset_float(file_id,datasetName, data);

    z.resize(numLayers+1);
    for (layer = 0; layer < numLayers+1; layer++) {
        z[layer].resize(numCol);
        for (col = 0; col < numCol; col++) {
            z[layer][col].resize(numRow,(int)0);
            for (row = 0; row < numRow; row++)
                z[layer][col][row] = -data[layer*(numCol)*(numRow)+col*(numRow)+row];
        }
    }

    // Load velocities
    getDatasetName("/Results/velocity U/velocity U_", index, datasetName);
    if (loadField3D(file_id, datasetName, u) != true) {
        return false;
    }
    getDatasetName("/Results/velocity V/velocity V_", index, datasetName);
    if (loadField3D(file_id, datasetName, v) != true) {
        return false;
    }
    getDatasetName("/Results/velocity W/velocity W_", index, datasetName);
    if (loadField3D(file_id, datasetName, w) != true) {
        return false;
    }

    // Load bathymetry and surface elevations as 2D arrays
    if (loadField2D(file_id, "/Grid/Bathymetry", bathymetry) != true) {
        return false;
    }
    getDatasetName("/Results/water level/water level_", index, datasetName);
    if (loadField2D(file_id, datasetName, surfaceElevation) != true) {
        return false;
    }

    // close file
    status = H5Fclose (file_id);

    delete [] data; data = NULL;
    delete [] maskData; maskData = NULL;
    return true;
}


//---------------------------------------------------------------------------------
//----------- loadAdditionalResults -----------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadAdditionalResults(char* hdf5FileName, int index)
{
    int col, row, layer;

    hid_t       file_id;            //Handle del archivo
    herr_t      status;             //Status de las funciones
    hsize_t     dims[3];            //Dimensiones del array
    char        datasetName[200];

    getDatasetName("", index, datasetName);
    cout << "Step " << datasetName << "- File " << hdf5FileName << " - Loading additional results..." << "\n";

    // Open HDF5 file
    file_id = H5Fopen (hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Get mask dimensions
    getDatasetName("/Grid/OpenPoints/OpenPoints_", index, datasetName);
    status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
    if (status == -1) return false;
    if (dims[0] != numLayers) return false;
    if (dims[1] != numCol) return false;
    if (dims[2] != numRow) return false;

    // Load additional fields
    getDatasetName("/Results/salinity/salinity_", index, datasetName);
    hasSalinity = loadField3D(file_id, datasetName, salinity);
    getDatasetName("/Results/density/density_", index, datasetName);
    hasDensity = loadField3D(file_id, datasetName, density);
    getDatasetName("/Results/fecal coliforms/fecal coliforms_", index, datasetName);
    hasColiforms = loadField3D(file_id, datasetName, coliforms);
    getDatasetName("/Results/temperature/temperature_", index, datasetName);
    hasTemperature = loadField3D(file_id, datasetName, temperature);

    //close file
    status = H5Fclose (file_id);

    return true;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

//---------------------------------------------------------------------------------
//----------- loadFieldsFromFieldFile ---------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadFieldsFromFieldFile(char* hdf5FileName, char* fieldsFileName, int index)
{
    int col, row, layer;

    hid_t       file_id;            //Handle del archivo
    herr_t      status;             //Status de las funciones
    hsize_t     dims[3];            //Dimensiones del array
    char        datasetName[200];

    getDatasetName("", index, datasetName);
    cout << "Step " << datasetName << "- File " << hdf5FileName << " - Loading requested fields..." << "\n";

    // Open HDF5 file
    file_id = H5Fopen (hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Get mask dimensions
    getDatasetName("/Grid/OpenPoints/OpenPoints_", index, datasetName);
    status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
    if (status == -1) return false;
    if (dims[0] != numLayers) return false;
    if (dims[1] != numCol) return false;
    if (dims[2] != numRow) return false;

    InputFile iFile(fieldsFileName);
    std::string line;
    while (iFile.returnNextLine(line)) {
        std::stringstream sstream;
        sstream << line << "_";
        getDatasetName(sstream.str().data(), index, datasetName);

        fields.resize(fields.size()+1);
        bool existe = loadField3D(file_id, datasetName, fields[fields.size()-1]);
        if (existe) {
            cout << "\tLoading field: " << line << "\n";
            std::vector<std::string> x = split(line, '/');
            fieldNames.push_back(x.back());
        } else {
            fields.resize(fields.size()-1);
        }
    }

    //close file
    status = H5Fclose (file_id);

    return true;
}

//---------------------------------------------------------------------------------
//----------- loadField3D ---------------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadField3D(hid_t    file_id, char* datasetName, vector<vector<vector<double> > > &outResultsVec)
{
    int col, row, layer;
    herr_t      status;         //Status de las funciones
    hsize_t     dims[3];        //Dimensiones del array
    float       *data;  //Datos

    data = (float *) malloc(sizeof(double) * (numRow+1) * (numCol+1) * (numLayers+1));

    // Get dataset dimensions
    status = H5LTget_dataset_info(file_id, datasetName,dims,NULL,NULL);
    if (dims[0] != numLayers) return false;
    if (dims[1] != numCol) return false;
    if (dims[2] != numRow) return false;

    // Read dataset
    status = H5LTread_dataset_float(file_id,datasetName, data);
    outResultsVec.resize(numLayers);
    for (layer = 0; layer < numLayers; layer++) {
        outResultsVec[layer].resize(numCol);
        for (col = 0; col < numCol; col++) {
            outResultsVec[layer][col].resize(numRow,(int)0);
            for (row = 0; row < numRow; row++) {
                outResultsVec[layer][col][row] = data[layer*(numCol)*(numRow)+col*(numRow)+row];
            }
        }
    }

    delete [] data; data = NULL;

    return true;
}

//---------------------------------------------------------------------------------
//----------- loadField2D ---------------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadField2D(hid_t    file_id, char* datasetName, vector<vector<double> > &outResultsVec)
{
    int col, row;
    herr_t      status;     //Status de las funciones
    hsize_t     dims[3];        //Dimensiones del array
    float       *data;      //Datos

    data = (float *) malloc(sizeof(double) * (numRow+1) * (numCol+1));

    // Get dataset dimensions
    status = H5LTget_dataset_info(file_id, datasetName,dims,NULL,NULL);
    if (dims[0] != numCol) return false;
    if (dims[1] != numRow) return false;

    // Read dataset
    status = H5LTread_dataset_float(file_id,datasetName, data);
    outResultsVec.resize(numCol);
    for (col = 0; col < numCol; col++) {
        outResultsVec[col].resize(numRow,(int)0);
        for (row = 0; row < numRow; row++) {
            outResultsVec[col][row] = data[col*(numRow)+row];
        }
    }

    delete [] data; data = NULL;

    return true;
}

//---------------------------------------------------------------------------------
//----------- loadMap -------------------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadMap( char* mapFileName )
{
    InputFile aFile( mapFileName );

    std::vector<std::string>    tagVec;

    aFile.setCommentPrefix(std::string("<"));

    map.resize(numCol);
    for (int col = 0; col < numCol; col++) {
        map[col].resize(numRow,0);
    }

    for (int row = 0; row < numRow; row++) {
        for (int col = 0; col < numCol; col++) {
            aFile.returnNextLineOfTags(tagVec);
            if (tagVec.size() < 3)
                return false;
            map[col][row] = strtod(tagVec[2].data(), 0);
        }
    }
    mapIsLoaded = true;
    return true;
}

//---------------------------------------------------------------------------------
//----------- convertResultsToVTK -------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::convertResultsToVTK(void)
{
    int row, col, layer;

    node.resize(0);
    cell.resize(0);
    numNodes = 0;
    numCells = 0;
    numValues = 0;

    cout << "\tConverting results to VTK..." << "\n";

    // Create array with nodes indexes
    cellLLBNodeIndex.resize(numLayers+1);
    for (layer = 0; layer < numLayers+1; layer++) {
        cellLLBNodeIndex[layer].resize(numCol+1);
        for (col=0; col<numCol+1; col++) {
            cellLLBNodeIndex[layer][col].resize(numRow+1, -1);
        }
    }

    for (layer=0; layer<numLayers; layer++) {
        for (col=0; col<numCol-4; col++) {
            for (row=0; row<numRow-4; row++) {
                if (mask[layer][col][row] == 1) {
                    // It's part of the domain

                    // Create cell and nodes if needed
                    cell.resize(cell.size()+1);
                    cell[numCells].node[0] = createNode(layer, col, row);
                    cell[numCells].node[1] = createNode(layer, col+1, row);
                    cell[numCells].node[2] = createNode(layer, col+1, row+1);
                    cell[numCells].node[3] = createNode(layer, col, row+1);
                    cell[numCells].node[4] = createNode(layer+1, col, row);
                    cell[numCells].node[5] = createNode(layer+1, col+1, row);
                    cell[numCells].node[6] = createNode(layer+1, col+1, row+1);
                    cell[numCells].node[7] = createNode(layer+1, col, row+1);
                    cell[numCells].numNodes = 8;
                    cell[numCells].col = col;
                    cell[numCells].row = row;
                    cell[numCells].layer = layer;
                    numCells++;
                    numValues += 9;
                }
            }
        }
    }
    return true;
}

//---------------------------------------------------------------------------------
//----------- createNode ----------------------------------------------------------
//---------------------------------------------------------------------------------
int MohidResults::createNode(int layer, int col, int row)
{
    int layerMask;
    if (layer < numLayers) {
        // It's a layer of cells
        layerMask = layer;
    } else {
        // It's the layer above the free surface
        layerMask = layer-1;
    }
    if (cellLLBNodeIndex[layer][col][row] == -1) {
        // Node doesn't exist --> create it
        double  sumaZ ;
        int     numSumaZ ;
        node.resize(node.size()+1);
        node[numNodes].x = x[col][row];
        node[numNodes].y = y[col][row];
        node[numNodes].layer = layer;
        sumaZ = 0;
        numSumaZ = 0;
        // Calculate average z coordinate of the neighbouring cells in rows and cols
        if (col > 0) {
            if (mask[layerMask][col-1][row] == 1) {
                numSumaZ++;
                sumaZ += z[layer][col-1][row];
            }
            if (row > 0) {
                if (mask[layerMask][col-1][row-1] == 1) {
                    numSumaZ++;
                    sumaZ += z[layer][col-1][row-1];
                }
            }
        }
        if (row > 0) {
            if (mask[layerMask][col][row-1] == 1) {
                numSumaZ++;
                sumaZ += z[layer][col][row-1];
            }
        }
        if (mask[layerMask][col][row] == 1) {
            numSumaZ++;
            sumaZ += z[layer][col][row];
        }

        node[numNodes].z = sumaZ / (double)numSumaZ;

        cellLLBNodeIndex[layer][col][row] = numNodes;
        numNodes++;
        return numNodes-1;
    } else {
        // Node exists
        return cellLLBNodeIndex[layer][col][row];
    }
}

//---------------------------------------------------------------------------------
//----------- writeResultsVTK -----------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::writeResultsVTK(char* vtkFileName, bool writeAs2D)
{
    vector<double>      nodeValues;
    vector<double>      nodeValues2;
    vector<double>      nodeValues3;
    nodeValues.resize(numNodes, 0);
    nodeValues2.resize(numNodes, 0);
    nodeValues3.resize(numNodes, 0);

    cout << "\tWriting VTK file: " << vtkFileName << "\n";

    // open file for output
    ofstream vtk(vtkFileName);

    // Header
    vtk << "# vtk DataFile Version 2.0" << "\n";
    vtk << vtkFileName << ", Created by Gmsh" << "\n";
    vtk << "ASCII" << "\n";
    vtk << "DATASET UNSTRUCTURED_GRID" << "\n";

    // nodes
    vtk << "POINTS " << numNodes << " double" << "\n";
    for (int i=0; i<numNodes; i++) {
        if (writeAs2D == true) {
            vtk << std::setprecision(10) << node[i].x << " " << node[i].y << " " << node[i].layer * 1.0 << "\n";
        } else {
            vtk << std::setprecision(10) << node[i].x << " " << node[i].y << " " << node[i].z  + verticalOffset << "\n";
        }
    }
    vtk << "\n";

    // cells
    vtk << "CELLS " << numCells << " " << numValues << "\n";
    for (int i=0; i<numCells; i++) {
        vtk << cell[i].numNodes << " " << cell[i].node[0]
                                    << " " << cell[i].node[1]
                                    << " " << cell[i].node[2]
                                    << " " << cell[i].node[3]
                                    << " " << cell[i].node[4]
                                    << " " << cell[i].node[5]
                                    << " " << cell[i].node[6]
                                    << " " << cell[i].node[7] << "\n";
    }
    vtk << "\n";
    vtk << "CELL_TYPES " << numCells << "\n";
    for (int i=0; i<numCells; i++) {
        vtk << "12" << "\n"; //VTK_HEXA
    }
    vtk << "\n";

    // Cell values  --------------------------------------
    vtk << "CELL_DATA " << numCells << "\n";
    // Velocity values
    vtk << "VECTORS " << "vel" << " float" << "\n";
    for (int i=0; i<numCells; i++) {
        vtk << u[cell[i].layer][cell[i].col][cell[i].row]
            << " " << v[cell[i].layer][cell[i].col][cell[i].row]
            << " " << w[cell[i].layer][cell[i].col][cell[i].row] << "\n";
    }
    vtk << "\n";

    //Valores de tirante
    vtk << "SCALARS " << "Tirante" << " float" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";
    for (int i=0; i<numCells; i++) {
        vtk << z[numLayers][cell[i].col][cell[i].row] - z[0][cell[i].col][cell[i].row] << "\n";
    }
    vtk << "\n";

    //Valores de fields adicionales
    for (int nc=0; nc<fieldNames.size(); nc++) {
        vtk << "SCALARS " << fieldNames[nc] << " float" << "\n";
        vtk << "LOOKUP_TABLE default" << "\n";
        for (int i=0; i<numCells; i++) {
            vtk << fields[nc][cell[i].layer][cell[i].col][cell[i].row] << "\n";
        }
        vtk << "\n";
    }

    // salinity
    if (hasSalinity) {
        vtk << "SCALARS " << "salinity" << " float" << "\n";
        vtk << "LOOKUP_TABLE default" << "\n";
        for (int i=0; i<numCells; i++) {
            vtk << salinity[cell[i].layer][cell[i].col][cell[i].row] << "\n";
        }
        vtk << "\n";
    }
    // density
    if (hasDensity) {
        vtk << "SCALARS " << "density" << " float" << "\n";
        vtk << "LOOKUP_TABLE default" << "\n";
        for (int i=0; i<numCells; i++) {
            vtk << density[cell[i].layer][cell[i].col][cell[i].row] << "\n";
        }
        vtk << "\n";
    }
    // coliforms
    if (hasColiforms) {
        vtk << "SCALARS " << "coliforms" << " float" << "\n";
        vtk << "LOOKUP_TABLE default" << "\n";
        for (int i=0; i<numCells; i++) {
            vtk << coliforms[cell[i].layer][cell[i].col][cell[i].row] << "\n";
        }
        vtk << "\n";
    }
    // temperature
    if (hasTemperature) {
        vtk << "SCALARS " << "temperature" << " float" << "\n";
        vtk << "LOOKUP_TABLE default" << "\n";
        for (int i=0; i<numCells; i++) {
            vtk << temperature[cell[i].layer][cell[i].col][cell[i].row] << "\n";
        }
        vtk << "\n";
    }
    // map values
    if (mapIsLoaded) {
        vtk << "SCALARS " << "map" << " float" << "\n";
        vtk << "LOOKUP_TABLE default" << "\n";
        for (int i=0; i<numCells; i++) {
            vtk << map[cell[i].col][cell[i].row] << "\n";
        }
        vtk << "\n";
    }

    // Node values --------------------------------------
    vtk << "POINT_DATA " << numNodes << "\n";
    // bathymetry
    calculateNodeValues2D(bathymetry, nodeValues);
    vtk << "SCALARS " << "bathymetry" << " float" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";
    for (int i=0; i<numNodes; i++) {
        vtk << verticalOffset - nodeValues[i] << "\n";
    }

    vtk << "\n";

    // surface elevation
    calculateNodeValues2D(surfaceElevation, nodeValues);
    vtk << "SCALARS " << "surfaceElevation" << " float" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";
    for (int i=0; i<numNodes; i++) {
        vtk << nodeValues[i] + verticalOffset << "\n";
    }

    vtk << "\n";

    // Velocity values
    calculateNodeValues(u, nodeValues);
    calculateNodeValues(v, nodeValues2);
    calculateNodeValues(w, nodeValues3);
    vtk << "VECTORS " << "vel" << " float" << "\n";
    for (int i=0; i<numNodes; i++) {
        vtk << nodeValues[i] << " " << nodeValues2[i] << " " << nodeValues3[i] << "\n";
    }

    vtk << "\n";
    vtk.close();


    return true;
}

//---------------------------------------------------------------------------------
//----------- writeResultsVTKBinary -----------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::writeResultsVTKBinary(char* vtkFileName, bool writeAs2D)
{
    int                 j,datoInt;
    float               dato;
    vector<double>      nodeValues;
    vector<double>      nodeValues2;
    vector<double>      nodeValues3;

    ifstream::pos_type  lastPos, newPos;

    nodeValues.resize(numNodes, 0);
    nodeValues2.resize(numNodes, 0);
    nodeValues3.resize(numNodes, 0);

    cout << "\tWriting VTK file: " << vtkFileName << "\n";

    // open file for output
    ofstream vtk(vtkFileName, ios::out | ios::binary);

    // Header
    vtk << "# vtk DataFile Version 2.0" << "\n";
    vtk << vtkFileName << ", Created by Gmsh" << "\n";
    vtk << "BINARY" << "\n";
    vtk << "DATASET UNSTRUCTURED_GRID" << "\n";

    cout << sizeof(int) << sizeof(float) << sizeof(double);
    // nodes
    vtk << "POINTS " << numNodes << " float" << "\n";

    lastPos = vtk.tellp();
    for (int i=0; i<numNodes; i++) {
        dato = node[i].x;
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
        dato = node[i].y;
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
        dato = node[i].z + verticalOffset;
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
        newPos = vtk.tellp();
        if (newPos-lastPos != 12) {
            cout << "\t\tERROR: Wrong node " << i << " - " << lastPos << " " << newPos << "\n";
            return false;
        }
        lastPos = newPos;
    }
    vtk << "\n";

    // cells
    vtk << "CELLS " << numCells << " " << numValues << "\n";
    lastPos = vtk.tellp();
    for (int i=0; i<numCells; i++) {
        datoInt = cell[i].numNodes;
        vtk.write(reinterpret_cast<const char *>(&datoInt)+3,sizeof(int)/4);
        vtk.write(reinterpret_cast<const char *>(&datoInt)+2,sizeof(int)/4);
        vtk.write(reinterpret_cast<const char *>(&datoInt)+1,sizeof(int)/4);
        vtk.write(reinterpret_cast<const char *>(&datoInt),sizeof(int)/4);

        for (int j=0; j<8; j++) {
            datoInt = cell[i].node[j];
            vtk.write(reinterpret_cast<const char *>(&datoInt)+3,sizeof(int)/4);
            vtk.write(reinterpret_cast<const char *>(&datoInt)+2,sizeof(int)/4);
            vtk.write(reinterpret_cast<const char *>(&datoInt)+1,sizeof(int)/4);
            vtk.write(reinterpret_cast<const char *>(&datoInt),sizeof(int)/4);
        }
        newPos = vtk.tellp();
        if (newPos-lastPos != 36) {
            cout << "\t\tERROR: Wrong cell " << i << " - " << lastPos << " " << newPos << "\n";
            return false;
        }
        lastPos = newPos;
    }
    vtk << "\n";
    vtk << "CELL_TYPES " << numCells << "\n";
    lastPos = vtk.tellp();
    for (int i=0; i<numCells; i++) {
        int a;
        datoInt = 12;
        vtk.write(reinterpret_cast<const char *>(&datoInt)+3,sizeof(int)/4);
        vtk.write(reinterpret_cast<const char *>(&datoInt)+2,sizeof(int)/4);
        vtk.write(reinterpret_cast<const char *>(&datoInt)+1,sizeof(int)/4);
        vtk.write(reinterpret_cast<const char *>(&datoInt),sizeof(int)/4);
        newPos = vtk.tellp();
        if (newPos-lastPos != 4) {
            cout << "\t\tERROR: Wrong cel_type " << i << " - " << lastPos << " " << newPos << "\n";
            return false;
        }
        lastPos = newPos;
        //vtk << "12" << "\n"; //VTK_HEXA
    }
    vtk << "\n";


    // Cell values  --------------------------------------
    vtk << "CELL_DATA " << numCells << "\n";
    // Velocity values
    vtk << "VECTORS " << "vel" << " float" << "\n";
    for (int i=0; i<numCells; i++) {
        //vtk << u[cell[i].layer][cell[i].col][cell[i].row]
        //  << " " << v[cell[i].layer][cell[i].col][cell[i].row]
        //  << " " << w[cell[i].layer][cell[i].col][cell[i].row] << "\n";
        dato = u[cell[i].layer][cell[i].col][cell[i].row];
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
        dato = v[cell[i].layer][cell[i].col][cell[i].row];
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
        dato = w[cell[i].layer][cell[i].col][cell[i].row];
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
    }
    vtk << "\n";

    // Node values --------------------------------------
    vtk << "POINT_DATA " << numNodes << "\n";
    // bathymetry
    calculateNodeValues2D(bathymetry, nodeValues);
    vtk << "SCALARS " << "bathymetry" << " float" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";
    for (int i=0; i<numNodes; i++) {
        dato = nodeValues[i] + verticalOffset;
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
    }
    vtk << "\n";
    vtk << "\n";

    // surface elevation
    calculateNodeValues2D(surfaceElevation, nodeValues);
    vtk << "SCALARS " << "surfaceElevation" << " float" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";
    for (int i=0; i<numNodes; i++) {
        dato = nodeValues[i] + verticalOffset;
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
    }
    vtk << "\n";
    vtk << "\n";

    // salinity
    if (hasSalinity) {
        calculateNodeValues(salinity, nodeValues);
        vtk << "SCALARS " << "salinity" << " float" << "\n";
        vtk << "LOOKUP_TABLE default" << "\n";
        for (int i=0; i<numNodes; i++) {
            dato = nodeValues[i];
            vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
            vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
            vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
            vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
        }
        vtk << "\n";
        vtk << "\n";
    }

    // Velocity values
    calculateNodeValues(u, nodeValues);
    calculateNodeValues(v, nodeValues2);
    calculateNodeValues(w, nodeValues3);
    vtk << "VECTORS " << "vel" << " float" << "\n";
    for (int i=0; i<numNodes; i++) {
        //vtk << nodeValues[i] << " " << nodeValues2[i] << " " << nodeValues3[i] << "\n";
        dato = nodeValues[i];
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
        dato = nodeValues2[i];
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
        dato = nodeValues3[i];
        vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
        vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
    }
    vtk << "\n";
    vtk << "\n";
    vtk.close();


    return true;
}


//---------------------------------------------------------------------------------
//----------- writeResultsASC -----------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::writeResultsASC(int index)
{
    char                fileName[500];

    vector<double>      nodeValues;
    vector<double>      nodeValues2;
    vector<double>      nodeValues3;

    nodeValues.resize(numNodes, 0);
    nodeValues2.resize(numNodes, 0);
    nodeValues3.resize(numNodes, 0);

    {
        getDatasetName("surfaceElev_", index, fileName);
        strcat ( fileName, ".asc");

        cout << "\tWriting VTK file: " << fileName << "\n";
        std::ofstream surfer(fileName);

        surfer << "ncols         " << numCol << "\n";
        surfer << "nrows         " << numRow << "\n";
        surfer << "xllcorner     " << x[0][0] << "\n";
        surfer << "yllcorner     " << y[0][0] << "\n";
        surfer << "cellsize      " << x[1][0] - x[0][0]<< "\n";
        surfer << "NODATA_value  -9999.0" << "\n";

        for (int k=numRow-1; k>=0; k--) {
            for (int i=0; i<numCol; i++) {
                if (mask[0][i][k] != 0) {
                    surfer << z[numLayers][i][k] + verticalOffset << " ";
                } else {
                    surfer << "-9999.0" << " ";
                }
            }
            surfer << "\n";
        }
        surfer.close();
    }
    {
        getDatasetName("bottomElev_", index, fileName);
        strcat ( fileName, ".asc");

        cout << "\tWriting VTK file: " << fileName << "\n";
        std::ofstream surfer(fileName);

        surfer << "ncols         " << numCol << "\n";
        surfer << "nrows         " << numRow << "\n";
        surfer << "xllcorner     " << x[0][0] << "\n";
        surfer << "yllcorner     " << y[0][0] << "\n";
        surfer << "cellsize      " << x[1][0] - x[0][0]<< "\n";
        surfer << "NODATA_value  -9999.0" << "\n";

        for (int k=numRow-1; k>=0; k--) {
            for (int i=0; i<numCol; i++) {
                if (mask[0][i][k] != 0) {
                    surfer << z[0][i][k] + verticalOffset << " ";
                } else {
                    surfer << "-9999.0" << " ";
                }
            }
            surfer << "\n";
        }
        surfer.close();
    }
    {
        getDatasetName("velociy_", index, fileName);
        strcat ( fileName, ".asc");

        cout << "\tWriting VTK file: " << fileName << "\n";
        std::ofstream surfer(fileName);

        surfer << "ncols         " << numCol << "\n";
        surfer << "nrows         " << numRow << "\n";
        surfer << "xllcorner     " << x[0][0] << "\n";
        surfer << "yllcorner     " << y[0][0] << "\n";
        surfer << "cellsize      " << x[1][0] - x[0][0]<< "\n";
        surfer << "NODATA_value  -9999.0" << "\n";

        for (int k=numRow-1; k>=0; k--) {
            for (int i=0; i<numCol; i++) {
                if (mask[0][i][k] != 0) {
                    double um = 0, vm = 0, wm = 0, modU = 0;
                    double tirante = 0;
                    for (int layer = 0; layer < numLayers; layer++) {
                        double dh = z[layer+1][i][k] - z[layer][i][k];
                        um += u[layer][i][k] * dh;
                        vm += v[layer][i][k] * dh;
                        wm += w[layer][i][k] * dh;
                        tirante += dh;
                    }
                    um /= tirante;
                    vm /= tirante;
                    wm /= tirante;
                    modU = sqrt(um*um+vm*vm+wm*wm);

                    surfer << modU << " ";
                } else {
                    surfer << "-9999.0" << " ";
                }
            }
            surfer << "\n";
        }
        surfer.close();
    }

    for (int nc=0; nc<fieldNames.size(); nc++) {
        getDatasetName(fieldNames[nc].data(), index, fileName);
        strcat ( fileName, ".asc");

        cout << "\tWriting VTK file: " << fileName << "\n";
        std::ofstream surfer(fileName);

        surfer << "ncols         " << numCol << "\n";
        surfer << "nrows         " << numRow << "\n";
        surfer << "xllcorner     " << x[0][0] << "\n";
        surfer << "yllcorner     " << y[0][0] << "\n";
        surfer << "cellsize      " << x[1][0] - x[0][0]<< "\n";
        surfer << "NODATA_value  -9999.0" << "\n";

        for (int k=numRow-1; k>=0; k--) {
            for (int i=0; i<numCol; i++) {
                if (mask[0][i][k] != 0) {
                    double medio = 0;
                    double tirante = 0;
                    for (int layer = 0; layer < numLayers; layer++) {
                        double dh = z[layer+1][i][k] - z[layer][i][k];
                        medio += fields[nc][layer][i][k] * dh;
                        tirante += dh;
                    }
                    medio /= tirante;
                    surfer << medio << " ";
                } else {
                    surfer << "-9999.0" << " ";
                }
            }
            surfer << "\n";
        }
        surfer.close();
    }


    return true;
}

//---------------------------------------------------------------------------------
//----------- calculateNodeValues -------------------------------------------------
//---------------------------------------------------------------------------------
void MohidResults::calculateNodeValues(vector<vector<vector<double> > > &inCellValues, vector<double> &outNodeValues)
{
    vector<int>         numNodeValues;

    numNodeValues.resize(outNodeValues.size(), 0);

    for (int i=0; i<numNodes; i++) {
        outNodeValues[i]=0;
    }

    for (int i=0; i<numCells; i++) {
        for (int j=0; j<cell[i].numNodes; j++) {
            outNodeValues[cell[i].node[j]] += inCellValues[cell[i].layer][cell[i].col][cell[i].row];
            numNodeValues[cell[i].node[j]]++;
        }
    }
    for (int i=0; i<numNodes; i++) {
        if (numNodeValues[i] > 0) {
            outNodeValues[i] = outNodeValues[i]/numNodeValues[i];
        } else {
            outNodeValues[i] = 0;
        }
    }
}

//---------------------------------------------------------------------------------
//----------- calculateNodeValues2D -----------------------------------------------
//---------------------------------------------------------------------------------
void MohidResults::calculateNodeValues2D(vector<vector<double> > &inCellValues, vector<double> &outNodeValues)
{
    vector<int>         numNodeValues;

    numNodeValues.resize(outNodeValues.size(), 0);

    for (int i=0; i<numNodes; i++) {
        outNodeValues[i]=0;
    }

    for (int i=0; i<numCells; i++) {
        for (int j=0; j<cell[i].numNodes; j++) {
            outNodeValues[cell[i].node[j]] += inCellValues[cell[i].col][cell[i].row];
            numNodeValues[cell[i].node[j]]++;
        }
    }
    for (int i=0; i<numNodes; i++) {
        if (numNodeValues[i] > 0) {
            outNodeValues[i] = outNodeValues[i]/numNodeValues[i];
        } else {
            outNodeValues[i] = 0;
        }
    }
}

//---------------------------------------------------------------------------------
//----------- getDatasetName ------------------------------------------------------
//---------------------------------------------------------------------------------
void MohidResults::getDatasetName(const char *baseName, int index, char *finalNameBfr)
{
    char            num[10];

    strcpy ( finalNameBfr, baseName );

    if (index < 10) {
        strcat ( finalNameBfr, "0000");
    } else if (index < 100){
        strcat ( finalNameBfr, "000");
    } else if (index < 1000){
        strcat ( finalNameBfr, "00");
    } else if (index < 10000){
        strcat ( finalNameBfr, "0");
    }

    sprintf(num, "%d", index);
    strcat ( finalNameBfr, num);
    return;
}
