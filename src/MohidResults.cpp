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
//----------- loadResult ----------------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadResult(char* archH5, int indice)
{
	int col, row, layer;

	hid_t       file_id;		//Handle del archivo
	herr_t      status;			//Status de las funciones
	hsize_t     dims[3];		//Dimensiones del array
	float       *data;			//Datos
	int			*mascara;//Mascara
	char		datasetName[200];

	cout << "\n";
	getDatasetName("", indice, datasetName);
	cout << "Archivo " << archH5 << " - Paso " << datasetName << " -- Cargando Resultados..." << "\n";

	// Abrir archivo
	file_id = H5Fopen (archH5, H5F_ACC_RDONLY, H5P_DEFAULT);

	// Obtener las dimensiones de la máscara
	getDatasetName("/Grid/OpenPoints/OpenPoints_", indice, datasetName);
	status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
	if (status == -1) return false;
	numLayers = dims[0];
	numCol = dims[1];
	numRow = dims[2];

	mascara = (int *)malloc(sizeof(int)*(numCol+1)*(numRow+1)*(numLayers+1));
	data = (float *)malloc(sizeof(float)*(numCol+1)*(numRow+1)*(numLayers+1));

	// Leer Dataset de máscara
	status = H5LTread_dataset_int(file_id,datasetName, mascara);

	mask.resize(numLayers);
	for (layer = 0; layer < numLayers; layer++) {
		mask[layer].resize(numCol);
		for (col = 0; col < numCol; col++) {
			mask[layer][col].resize(numRow,(int)0);
			for (row = 0; row < numRow; row++) {
				mask[layer][col][row] = mascara[layer*numCol*numRow+col*numRow+row];
			}
		}
	}


	// Obtener las dimensiones de ConnectionX
	status = H5LTget_dataset_info(file_id,"/Grid/ConnectionX",dims,NULL,NULL);
	if (status == -1) return false;
	if (dims[0] != numCol+1) return false;
	if (dims[1] != numRow+1) return false;

	// Leer Dataset de ConnectionX
	status = H5LTread_dataset_float(file_id,"/Grid/ConnectionX", data);

	x.resize(numCol+1);
	for (col = 0; col < numCol+1; col++) {
		x[col].resize(numRow+1,(float)0);
		for (row = 0; row < numRow+1; row++) {
			x[col][row] = data[col*(numRow+1)+row];
		}
	}


	// Obtener las dimensiones de ConnectionY
	status = H5LTget_dataset_info(file_id,"/Grid/ConnectionY",dims,NULL,NULL);
	if (status == -1) return false;
	if (dims[0] != numCol+1) return false;
	if (dims[1] != numRow+1) return false;

	// Leer Dataset de ConnectionY
	status = H5LTread_dataset_float(file_id,"/Grid/ConnectionY", data);

	y.resize(numCol+1);
	for (col = 0; col < numCol+1; col++) {
		y[col].resize(numRow+1,(float)0);
		for (row = 0; row < numRow+1; row++) {
			y[col][row] = data[col*(numRow+1)+row];
		}
	}


	// Obtener las dimensiones de VerticalZ
	getDatasetName("/Grid/VerticalZ/Vertical_", indice, datasetName);
	status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
	if (status == -1) return false;
	if (dims[0] != numLayers+1) return false;
	if (dims[1] != numCol) return false;
	if (dims[2] != numRow) return false;

	// Leer Dataset de VerticalZ
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

	// Cargar Velocidades
	getDatasetName("/Results/velocity U/velocity U_", indice, datasetName);
	if (loadField3D(file_id, datasetName, u) != true) {
		return false;
	}
	getDatasetName("/Results/velocity V/velocity V_", indice, datasetName);
	if (loadField3D(file_id, datasetName, v) != true) {
		return false;
	}
	getDatasetName("/Results/velocity W/velocity W_", indice, datasetName);
	if (loadField3D(file_id, datasetName, w) != true) {
		return false;
	}

	// Cargar bathymetry, y nivel superficial cmo arrays planos (columna, fila)
	if (loadField2D(file_id, "/Grid/Bathymetry", bathymetry) != true) {
		return false;
	}

	getDatasetName("/Results/water level/water level_", indice, datasetName);
	if (loadField2D(file_id, datasetName, surfaceElevation) != true) {
		return false;
	}

	//close file
	status = H5Fclose (file_id);

	delete [] data; data = NULL;
	delete [] mascara; mascara = NULL;
	return true;
}


//---------------------------------------------------------------------------------
//----------- addResult -----------------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::addResult(char* archH5, int indice)
{
	int col, row, layer;

	hid_t       file_id;			//Handle del archivo
	herr_t      status;				//Status de las funciones
	hsize_t     dims[3];			//Dimensiones del array
	char		datasetName[200];

	getDatasetName("", indice, datasetName);
	cout << "Archivo " << archH5 << " - Paso " << datasetName << " -- Cargando Resultados..." << "\n";

	// Abrir archivo
	file_id = H5Fopen (archH5, H5F_ACC_RDONLY, H5P_DEFAULT);

	// Obtener las dimensiones de la máscara
	getDatasetName("/Grid/OpenPoints/OpenPoints_", indice, datasetName);
	status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
	if (status == -1) return false;
	if (dims[0] != numLayers) return false;
	if (dims[1] != numCol) return false;
	if (dims[2] != numRow) return false;

	// Cargar Campos de calidad
	getDatasetName("/Results/salinity/salinity_", indice, datasetName);
	hasSalinity = loadField3D(file_id, datasetName, salinity);

	getDatasetName("/Results/density/density_", indice, datasetName);
	hasDensity = loadField3D(file_id, datasetName, density);

	getDatasetName("/Results/fecal coliforms/fecal coliforms_", indice, datasetName);
	hasColiforms = loadField3D(file_id, datasetName, coliforms);

	getDatasetName("/Results/temperature/temperature_", indice, datasetName);
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
bool MohidResults::loadFieldsFromFieldFile(char* archH5, char* archCampos, int indice)
{
	int col, row, layer;

	hid_t       file_id;			//Handle del archivo
	herr_t      status;				//Status de las funciones
	hsize_t     dims[3];			//Dimensiones del array
	char		datasetName[200];

	getDatasetName("", indice, datasetName);
	cout << "Archivo " << archH5 << " - Paso " << datasetName << " -- Cargando Campos Adicionales..." << "\n";

	// Abrir archivo
	file_id = H5Fopen (archH5, H5F_ACC_RDONLY, H5P_DEFAULT);

	// Obtener las dimensiones de la máscara
	getDatasetName("/Grid/OpenPoints/OpenPoints_", indice, datasetName);
	status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
	if (status == -1) return false;
	if (dims[0] != numLayers) return false;
	if (dims[1] != numCol) return false;
	if (dims[2] != numRow) return false;

	InputFile iFile(archCampos);
	std::string line;
	while (iFile.returnNextLine(line)) {
		std::stringstream sstream;
		sstream << line << "_";
		getDatasetName(sstream.str().data(), indice, datasetName);

		fields.resize(fields.size()+1);
		bool existe = loadField3D(file_id, datasetName, fields[fields.size()-1]);
		if (existe) {
			cout << "\tCargado campo: " << line << "\n";

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
bool MohidResults::loadField3D(hid_t	file_id, char* nombreDataset, vector<vector<vector<double> > > &vectorRes)
{
	int col, row, layer;
	herr_t      status;			//Status de las funciones
	hsize_t     dims[3];		//Dimensiones del array
	float       *data;	//Datos

	data = (float *) malloc(sizeof(double) * (numRow+1) * (numCol+1) * (numLayers+1));

	// Obtener las dimensiones del dataset
	status = H5LTget_dataset_info(file_id, nombreDataset,dims,NULL,NULL);
	if (dims[0] != numLayers) return false;
	if (dims[1] != numCol) return false;
	if (dims[2] != numRow) return false;

	// Leer Dataset
	status = H5LTread_dataset_float(file_id,nombreDataset, data);

	vectorRes.resize(numLayers);
	for (layer = 0; layer < numLayers; layer++) {
		vectorRes[layer].resize(numCol);
		for (col = 0; col < numCol; col++) {
			vectorRes[layer][col].resize(numRow,(int)0);
			for (row = 0; row < numRow; row++) {
				vectorRes[layer][col][row] = data[layer*(numCol)*(numRow)+col*(numRow)+row];
			}
		}
	}

	delete [] data; data = NULL;

	return true;
}

//---------------------------------------------------------------------------------
//----------- loadField2D ---------------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadField2D(hid_t	file_id, char* nombreDataset, vector<vector<double> > &vectorRes)
{
	int col, row;
	herr_t      status;		//Status de las funciones
	hsize_t     dims[3];		//Dimensiones del array
	float       *data;		//Datos

	data = (float *) malloc(sizeof(double) * (numRow+1) * (numCol+1));

	// Obtener las dimensiones del dataset
	status = H5LTget_dataset_info(file_id, nombreDataset,dims,NULL,NULL);
	if (dims[0] != numCol) return false;
	if (dims[1] != numRow) return false;

	// Leer Dataset
	status = H5LTread_dataset_float(file_id,nombreDataset, data);

	vectorRes.resize(numCol);
	for (col = 0; col < numCol; col++) {
		vectorRes[col].resize(numRow,(int)0);
		for (row = 0; row < numRow; row++) {
			vectorRes[col][row] = data[col*(numRow)+row];
		}
	}

	delete [] data; data = NULL;

	return true;
}

//---------------------------------------------------------------------------------
//----------- loadMap -------------------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::loadMap( char* nombreMapa )
{
	InputFile aFile( nombreMapa );

	std::vector<std::string>	tagVec;

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
	int	row, col, layer;

	node.resize(0);
	cell.resize(0);
	numNodes = 0;
	numCells = 0;
	numValues = 0;

	cout << "    Convirtiendo los resultados a vtk..." << "\n";

	//Crear array de número de nodes
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
					//Es parte del dominio

					//Crear el cell y los Nodos si hace falta
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
		// Es una layer de celdas
		layerMask = layer;
	} else {
		// Es la layer que está por encima de la superficie libre
		layerMask = layer-1;
	}
	if (cellLLBNodeIndex[layer][col][row] == -1) {
		//El node no existe --> crearlo
		double	sumaZ ;
		int		numSumaZ ;
		node.resize(node.size()+1);
		node[numNodes].x = x[col][row];
		node[numNodes].y = y[col][row];
		node[numNodes].layer = layer;
		sumaZ = 0;
		numSumaZ = 0;
		//Calcular la coordenada Z promedio de las celdas vecinas en fias y columnas
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
		//El node existe
		return cellLLBNodeIndex[layer][col][row];
	}
}

//---------------------------------------------------------------------------------
//----------- writeResultsVTK -----------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::writeResultsVTK(char* archVTK, bool EscribirComo2D)
{
	int					i;
	vector<double>		valorNodo;
	vector<double>		valorNodo2;
	vector<double>		valorNodo3;

	valorNodo.resize(numNodes, 0);
	valorNodo2.resize(numNodes, 0);
	valorNodo3.resize(numNodes, 0);

	cout << "    Escribiendo los resultados en el archivo vtk..." << "\n";

	// open file for output
	ofstream vtk(archVTK);

	//Encabezado
	vtk << "# vtk DataFile Version 2.0" << "\n";
	vtk << archVTK << ", Created by Gmsh" << "\n";
	vtk << "ASCII" << "\n";
	vtk << "DATASET UNSTRUCTURED_GRID" << "\n";

	//Nodos
	vtk << "POINTS " << numNodes << " double" << "\n";
	for (i=0; i<numNodes; i++) {
		if (EscribirComo2D == true) {
			vtk << std::setprecision(10) << node[i].x << " " << node[i].y << " " << node[i].layer * 1.0 << "\n";
		} else {
			vtk << std::setprecision(10) << node[i].x << " " << node[i].y << " " << node[i].z  + verticalOffset << "\n";
		}
	}
	vtk << "\n";

	//Celdas
	vtk << "CELLS " << numCells << " " << numValues << "\n";
	for (i=0; i<numCells; i++) {
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
	for (i=0; i<numCells; i++) {
		vtk << "12" << "\n"; //VTK_HEXA
	}
	vtk << "\n";

	//Valores en Celdas --------------------------------------
	vtk << "CELL_DATA " << numCells << "\n";
	//Valores de Velocidad
	vtk << "VECTORS " << "vel" << " float" << "\n";
	for (i=0; i<numCells; i++) {
		vtk << u[cell[i].layer][cell[i].col][cell[i].row]
			<< " " << v[cell[i].layer][cell[i].col][cell[i].row]
			<< " " << w[cell[i].layer][cell[i].col][cell[i].row] << "\n";
	}
	vtk << "\n";

	//Valores de tirante
	vtk << "SCALARS " << "Tirante" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numCells; i++) {
		vtk << z[numLayers][cell[i].col][cell[i].row] - z[0][cell[i].col][cell[i].row] << "\n";
	}
	vtk << "\n";

	//Valores de fields adicionales
	for (int nc=0; nc<fieldNames.size(); nc++) {
		vtk << "SCALARS " << fieldNames[nc] << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numCells; i++) {
			vtk << fields[nc][cell[i].layer][cell[i].col][cell[i].row] << "\n";
		}
		vtk << "\n";
	}

	//Valores de salinidad
	if (hasSalinity) {
		vtk << "SCALARS " << "salinity" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numCells; i++) {
			vtk << salinity[cell[i].layer][cell[i].col][cell[i].row] << "\n";
		}
		vtk << "\n";
	}
	//Valores de densidad
	if (hasDensity) {
		vtk << "SCALARS " << "density" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numCells; i++) {
			vtk << density[cell[i].layer][cell[i].col][cell[i].row] << "\n";
		}
		vtk << "\n";
	}
	//Valores de coliformes
	if (hasColiforms) {
		vtk << "SCALARS " << "coliforms" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numCells; i++) {
			vtk << coliforms[cell[i].layer][cell[i].col][cell[i].row] << "\n";
		}
		vtk << "\n";
	}
	//Valores de temperatura
	if (hasTemperature) {
		vtk << "SCALARS " << "temperature" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numCells; i++) {
			vtk << temperature[cell[i].layer][cell[i].col][cell[i].row] << "\n";
		}
		vtk << "\n";
	}
	//Valores del map
	if (mapIsLoaded) {
		vtk << "SCALARS " << "Mapa" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numCells; i++) {
			vtk << map[cell[i].col][cell[i].row] << "\n";
		}
		vtk << "\n";
	}

	// Valores en Nodos --------------------------------------
	vtk << "POINT_DATA " << numNodes << "\n";
	//Valores de bathymetry
	calculateNodeValues2D(bathymetry, valorNodo);
	vtk << "SCALARS " << "bathymetry" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numNodes; i++) {
		vtk << verticalOffset - valorNodo[i] << "\n";
	}

	vtk << "\n";

	//Valores de nivel de agua
	calculateNodeValues2D(surfaceElevation, valorNodo);
	vtk << "SCALARS " << "surfaceElevation" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numNodes; i++) {
		vtk << valorNodo[i] + verticalOffset << "\n";
	}

	vtk << "\n";

	//Valores de Velocidad
	calculateNodeValues(u, valorNodo);
	calculateNodeValues(v, valorNodo2);
	calculateNodeValues(w, valorNodo3);
	vtk << "VECTORS " << "vel" << " float" << "\n";
	for (i=0; i<numNodes; i++) {
		vtk << valorNodo[i] << " " << valorNodo2[i] << " " << valorNodo3[i] << "\n";
	}

	vtk << "\n";
	vtk.close();

	cout << "    Listo!" << "\n";

	return true;
}

//---------------------------------------------------------------------------------
//----------- writeResultsVTKBinary -----------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::writeResultsVTKBinary(char* archVTK, bool EscribirComo2D)
{
	int					i,j,datoInt;
	float				dato;
	vector<double>		valorNodo;
	vector<double>		valorNodo2;
	vector<double>		valorNodo3;

	ifstream::pos_type	lastPos, newPos;

	valorNodo.resize(numNodes, 0);
	valorNodo2.resize(numNodes, 0);
	valorNodo3.resize(numNodes, 0);

	cout << "    Escribiendo los resultados en el archivo vtk..." << "\n";

	// open file for output
	ofstream vtk(archVTK, ios::out | ios::binary);

	//Encabezado
	vtk << "# vtk DataFile Version 2.0" << "\n";
	vtk << archVTK << ", Created by Gmsh" << "\n";
	vtk << "BINARY" << "\n";
	vtk << "DATASET UNSTRUCTURED_GRID" << "\n";

	cout <<	sizeof(int) << sizeof(float) << sizeof(double);
	//Nodos
	vtk << "POINTS " << numNodes << " float" << "\n";

	lastPos = vtk.tellp();
	for (i=0; i<numNodes; i++) {
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
			cout << "    Mal node " << i << " - " << lastPos << " " << newPos << "\n";
			return false;
		}
		lastPos = newPos;
	}
	vtk << "\n";
	newPos = vtk.tellp();
	if (newPos-lastPos != 1) {
		cout << "    Mal enter 1 " << " - " << lastPos << " " << newPos << "\n";
		return false;
	}
	lastPos = newPos;

	//Celdas
	vtk << "CELLS " << numCells << " " << numValues << "\n";
	lastPos = vtk.tellp();
	for (i=0; i<numCells; i++) {
		datoInt = cell[i].numNodes;
		vtk.write(reinterpret_cast<const char *>(&datoInt)+3,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt)+2,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt)+1,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt),sizeof(int)/4);

		for (j=0; j<8; j++) {
			datoInt = cell[i].node[j];
			vtk.write(reinterpret_cast<const char *>(&datoInt)+3,sizeof(int)/4);
			vtk.write(reinterpret_cast<const char *>(&datoInt)+2,sizeof(int)/4);
			vtk.write(reinterpret_cast<const char *>(&datoInt)+1,sizeof(int)/4);
			vtk.write(reinterpret_cast<const char *>(&datoInt),sizeof(int)/4);
		}
		newPos = vtk.tellp();
		if (newPos-lastPos != 36) {
			cout << "    Mal cell " << i << " - " << lastPos << " " << newPos << "\n";
			return false;
		}
		lastPos = newPos;
	}
	vtk << "\n";
	vtk << "CELL_TYPES " << numCells << "\n";
	lastPos = vtk.tellp();
	for (i=0; i<numCells; i++) {
		int a;
		datoInt = 12;
		vtk.write(reinterpret_cast<const char *>(&datoInt)+3,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt)+2,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt)+1,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt),sizeof(int)/4);
		newPos = vtk.tellp();
		if (newPos-lastPos != 4) {
			cout << "    Mal cel_type " << i << " - " << lastPos << " " << newPos << "\n";
			return false;
		}
		lastPos = newPos;
		//vtk << "12" << "\n"; //VTK_HEXA
	}
	vtk << "\n";


	//Valores en Celdas --------------------------------------
	vtk << "CELL_DATA " << numCells << "\n";
	//Valores de Velocidad
	vtk << "VECTORS " << "vel" << " float" << "\n";
	for (i=0; i<numCells; i++) {
		//vtk << u[cell[i].layer][cell[i].col][cell[i].row]
		//	<< " " << v[cell[i].layer][cell[i].col][cell[i].row]
		//	<< " " << w[cell[i].layer][cell[i].col][cell[i].row] << "\n";
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

	// Valores en Nodos --------------------------------------
	vtk << "POINT_DATA " << numNodes << "\n";
	//Valores de bathymetry
	calculateNodeValues2D(bathymetry, valorNodo);
	vtk << "SCALARS " << "bathymetry" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numNodes; i++) {
		dato = valorNodo[i] + verticalOffset;
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
	}
	vtk << "\n";
	vtk << "\n";

	//Valores de nivel de agua
	calculateNodeValues2D(surfaceElevation, valorNodo);
	vtk << "SCALARS " << "surfaceElevation" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numNodes; i++) {
		dato = valorNodo[i] + verticalOffset;
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
	}
	vtk << "\n";
	vtk << "\n";

	//Valores de salinidad
	if (hasSalinity) {
		calculateNodeValues(salinity, valorNodo);
		vtk << "SCALARS " << "salinity" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numNodes; i++) {
			dato = valorNodo[i];
			vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
			vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
			vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
			vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
		}
		vtk << "\n";
		vtk << "\n";
	}

	//Valores de Velocidad
	calculateNodeValues(u, valorNodo);
	calculateNodeValues(v, valorNodo2);
	calculateNodeValues(w, valorNodo3);
	vtk << "VECTORS " << "vel" << " float" << "\n";
	for (i=0; i<numNodes; i++) {
		//vtk << valorNodo[i] << " " << valorNodo2[i] << " " << valorNodo3[i] << "\n";
		dato = valorNodo[i];
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
		dato = valorNodo2[i];
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
		dato = valorNodo3[i];
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
	}
	vtk << "\n";
	vtk << "\n";
	vtk.close();

	cout << "    Listo!" << "\n";

	return true;
}


//---------------------------------------------------------------------------------
//----------- writeResultsASC -----------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::writeResultsASC(int indice)
{
	char				nombreArchivo[500];

	vector<double>		valorNodo;
	vector<double>		valorNodo2;
	vector<double>		valorNodo3;

	valorNodo.resize(numNodes, 0);
	valorNodo2.resize(numNodes, 0);
	valorNodo3.resize(numNodes, 0);

	cout << "    Escribiendo los resultados en formato del GIS..." << "\n";
	{
		getDatasetName("supLibre_", indice, nombreArchivo);
		strcat ( nombreArchivo, ".asc");

		cout << "\t\tEscribiendo " << nombreArchivo << "\n";
		std::ofstream surfer(nombreArchivo);

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
		getDatasetName("nivelFondo_", indice, nombreArchivo);
		strcat ( nombreArchivo, ".asc");

		cout << "\t\tEscribiendo " << nombreArchivo << "\n";
		std::ofstream surfer(nombreArchivo);

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
		getDatasetName("moduloVel_", indice, nombreArchivo);
		strcat ( nombreArchivo, ".asc");

		cout << "\t\tEscribiendo " << nombreArchivo << "\n";
		std::ofstream surfer(nombreArchivo);

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
		getDatasetName(fieldNames[nc].data(), indice, nombreArchivo);
		strcat ( nombreArchivo, ".asc");

		cout << "\t\tEscribiendo " << nombreArchivo << "\n";
		std::ofstream surfer(nombreArchivo);

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

	cout << "\tListo!" << "\n";

	return true;
}

//---------------------------------------------------------------------------------
//----------- calculateNodeValues -------------------------------------------------
//---------------------------------------------------------------------------------
void MohidResults::calculateNodeValues(vector<vector<vector<double> > > &valElemento, vector<double> &valNodo)
{
	int					i, j;
	vector<int>			numValNodo;

	numValNodo.resize(valNodo.size(), 0);

	for (i=0; i<numNodes; i++) {
		valNodo[i]=0;
	}

	for (i=0; i<numCells; i++) {
		for (j=0; j<cell[i].numNodes; j++) {
			valNodo[cell[i].node[j]] += valElemento[cell[i].layer][cell[i].col][cell[i].row];
			numValNodo[cell[i].node[j]]++;
		}
	}
	for (i=0; i<numNodes; i++) {
		if (numValNodo[i] > 0) {
			valNodo[i] = valNodo[i]/numValNodo[i];
		} else {
			valNodo[i] = 0;
		}
	}
}

//---------------------------------------------------------------------------------
//----------- calculateNodeValues2D -----------------------------------------------
//---------------------------------------------------------------------------------
void MohidResults::calculateNodeValues2D(vector<vector<double> > &valElemento, vector<double> &valNodo)
{
	int					i, j;
	vector<int>			numValNodo;

	numValNodo.resize(valNodo.size(), 0);

	for (i=0; i<numNodes; i++) {
		valNodo[i]=0;
	}

	for (i=0; i<numCells; i++) {
		for (j=0; j<cell[i].numNodes; j++) {
			valNodo[cell[i].node[j]] += valElemento[cell[i].col][cell[i].row];
			numValNodo[cell[i].node[j]]++;
		}
	}
	for (i=0; i<numNodes; i++) {
		if (numValNodo[i] > 0) {
			valNodo[i] = valNodo[i]/numValNodo[i];
		} else {
			valNodo[i] = 0;
		}
	}
}

//---------------------------------------------------------------------------------
//----------- getDatasetName ------------------------------------------------------
//---------------------------------------------------------------------------------
void MohidResults::getDatasetName(const char *nombreBase, int indice, char *nombreFinal)
{
	char			num[10];

	strcpy ( nombreFinal, nombreBase );

	if (indice < 10) {
		strcat ( nombreFinal, "0000");
	} else if (indice < 100){
		strcat ( nombreFinal, "000");
	} else if (indice < 1000){
		strcat ( nombreFinal, "00");
	} else if (indice < 10000){
		strcat ( nombreFinal, "0");
	}

	sprintf(num, "%d", indice);
	strcat ( nombreFinal, num);
	return;
}
