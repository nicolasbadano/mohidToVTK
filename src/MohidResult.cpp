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
	offsetVertical = offset;
	cargadaMalla = false;
	cargadoResultado = false;
	cargadoMapa = false;

	campos.reserve(5);

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
	int col, fil, capa;

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
	numCapas = dims[0];
	numCol = dims[1];
	numFil = dims[2];

	mascara = (int *)malloc(sizeof(int)*(numCol+1)*(numFil+1)*(numCapas+1));
	data = (float *)malloc(sizeof(float)*(numCol+1)*(numFil+1)*(numCapas+1));

	// Leer Dataset de máscara
	status = H5LTread_dataset_int(file_id,datasetName, mascara);

	mask.resize(numCapas);
	for (capa = 0; capa < numCapas; capa++) {
		mask[capa].resize(numCol);
		for (col = 0; col < numCol; col++) {
			mask[capa][col].resize(numFil,(int)0);
			for (fil = 0; fil < numFil; fil++) {
				mask[capa][col][fil] = mascara[capa*numCol*numFil+col*numFil+fil];
			}
		}
	}


	// Obtener las dimensiones de ConnectionX
	status = H5LTget_dataset_info(file_id,"/Grid/ConnectionX",dims,NULL,NULL);
	if (status == -1) return false;
	if (dims[0] != numCol+1) return false;
	if (dims[1] != numFil+1) return false;

	// Leer Dataset de ConnectionX
	status = H5LTread_dataset_float(file_id,"/Grid/ConnectionX", data);

	x.resize(numCol+1);
	for (col = 0; col < numCol+1; col++) {
		x[col].resize(numFil+1,(float)0);
		for (fil = 0; fil < numFil+1; fil++) {
			x[col][fil] = data[col*(numFil+1)+fil];
		}
	}


	// Obtener las dimensiones de ConnectionY
	status = H5LTget_dataset_info(file_id,"/Grid/ConnectionY",dims,NULL,NULL);
	if (status == -1) return false;
	if (dims[0] != numCol+1) return false;
	if (dims[1] != numFil+1) return false;

	// Leer Dataset de ConnectionY
	status = H5LTread_dataset_float(file_id,"/Grid/ConnectionY", data);

	y.resize(numCol+1);
	for (col = 0; col < numCol+1; col++) {
		y[col].resize(numFil+1,(float)0);
		for (fil = 0; fil < numFil+1; fil++) {
			y[col][fil] = data[col*(numFil+1)+fil];
		}
	}


	// Obtener las dimensiones de VerticalZ
	getDatasetName("/Grid/VerticalZ/Vertical_", indice, datasetName);
	status = H5LTget_dataset_info(file_id,datasetName,dims,NULL,NULL);
	if (status == -1) return false;
	if (dims[0] != numCapas+1) return false;
	if (dims[1] != numCol) return false;
	if (dims[2] != numFil) return false;

	// Leer Dataset de VerticalZ
	status = H5LTread_dataset_float(file_id,datasetName, data);

	z.resize(numCapas+1);
	for (capa = 0; capa < numCapas+1; capa++) {
		z[capa].resize(numCol);
		for (col = 0; col < numCol; col++) {
			z[capa][col].resize(numFil,(int)0);
			for (fil = 0; fil < numFil; fil++)
				z[capa][col][fil] = -data[capa*(numCol)*(numFil)+col*(numFil)+fil];
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

	// Cargar batimetria, y nivel superficial cmo arrays planos (columna, fila)
	if (loadField2D(file_id, "/Grid/Bathymetry", batimetria) != true) {
		return false;
	}

	getDatasetName("/Results/water level/water level_", indice, datasetName);
	if (loadField2D(file_id, datasetName, nivelSuperficie) != true) {
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
	int col, fil, capa;

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
	if (dims[0] != numCapas) return false;
	if (dims[1] != numCol) return false;
	if (dims[2] != numFil) return false;

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
	int col, fil, capa;

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
	if (dims[0] != numCapas) return false;
	if (dims[1] != numCol) return false;
	if (dims[2] != numFil) return false;

	InputFile iFile(archCampos);
	std::string line;
	while (iFile.returnNextLine(line)) {
		std::stringstream sstream;
		sstream << line << "_";
		getDatasetName(sstream.str().data(), indice, datasetName);

		campos.resize(campos.size()+1);
		bool existe = loadField3D(file_id, datasetName, campos[campos.size()-1]);
		if (existe) {
			cout << "\tCargado campo: " << line << "\n";

			std::vector<std::string> x = split(line, '/');

			listaCampos.push_back(x.back());
		} else {
			campos.resize(campos.size()-1);
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
	int col, fil, capa;
	herr_t      status;			//Status de las funciones
	hsize_t     dims[3];		//Dimensiones del array
	float       *data;	//Datos

	data = (float *) malloc(sizeof(double) * (numFil+1) * (numCol+1) * (numCapas+1));

	// Obtener las dimensiones del dataset
	status = H5LTget_dataset_info(file_id, nombreDataset,dims,NULL,NULL);
	if (dims[0] != numCapas) return false;
	if (dims[1] != numCol) return false;
	if (dims[2] != numFil) return false;

	// Leer Dataset
	status = H5LTread_dataset_float(file_id,nombreDataset, data);

	vectorRes.resize(numCapas);
	for (capa = 0; capa < numCapas; capa++) {
		vectorRes[capa].resize(numCol);
		for (col = 0; col < numCol; col++) {
			vectorRes[capa][col].resize(numFil,(int)0);
			for (fil = 0; fil < numFil; fil++) {
				vectorRes[capa][col][fil] = data[capa*(numCol)*(numFil)+col*(numFil)+fil];
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
	int col, fil;
	herr_t      status;		//Status de las funciones
	hsize_t     dims[3];		//Dimensiones del array
	float       *data;		//Datos

	data = (float *) malloc(sizeof(double) * (numFil+1) * (numCol+1));

	// Obtener las dimensiones del dataset
	status = H5LTget_dataset_info(file_id, nombreDataset,dims,NULL,NULL);
	if (dims[0] != numCol) return false;
	if (dims[1] != numFil) return false;

	// Leer Dataset
	status = H5LTread_dataset_float(file_id,nombreDataset, data);

	vectorRes.resize(numCol);
	for (col = 0; col < numCol; col++) {
		vectorRes[col].resize(numFil,(int)0);
		for (fil = 0; fil < numFil; fil++) {
			vectorRes[col][fil] = data[col*(numFil)+fil];
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

	mapa.resize(numCol);
	for (int col = 0; col < numCol; col++) {
		mapa[col].resize(numFil,0);
	}

	for (int fil = 0; fil < numFil; fil++) {
		for (int col = 0; col < numCol; col++) {
			aFile.returnNextLineOfTags(tagVec);
			if (tagVec.size() < 3)
				return false;
			mapa[col][fil] = strtod(tagVec[2].data(), 0);
		}
	}
	cargadoMapa = true;
	return true;
}

//---------------------------------------------------------------------------------
//----------- convertResultsToVTK -------------------------------------------------
//---------------------------------------------------------------------------------
bool MohidResults::convertResultsToVTK(void)
{
	int	fil, col, capa;

	nodo.resize(0);
	elemento.resize(0);
	numNodos = 0;
	numElementos = 0;
	numDatos = 0;

	cout << "    Convirtiendo los resultados a vtk..." << "\n";

	//Crear array de número de nodos
	nodollbCelda.resize(numCapas+1);
	for (capa = 0; capa < numCapas+1; capa++) {
		nodollbCelda[capa].resize(numCol+1);
		for (col=0; col<numCol+1; col++) {
			nodollbCelda[capa][col].resize(numFil+1, -1);
		}
	}

	for (capa=0; capa<numCapas; capa++) {
		for (col=0; col<numCol-4; col++) {
			for (fil=0; fil<numFil-4; fil++) {
				if (mask[capa][col][fil] == 1) {
					//Es parte del dominio

					//Crear el elemento y los Nodos si hace falta
					elemento.resize(elemento.size()+1);
					elemento[numElementos].nodo[0] = createNode(capa, col, fil);
					elemento[numElementos].nodo[1] = createNode(capa, col+1, fil);
					elemento[numElementos].nodo[2] = createNode(capa, col+1, fil+1);
					elemento[numElementos].nodo[3] = createNode(capa, col, fil+1);
					elemento[numElementos].nodo[4] = createNode(capa+1, col, fil);
					elemento[numElementos].nodo[5] = createNode(capa+1, col+1, fil);
					elemento[numElementos].nodo[6] = createNode(capa+1, col+1, fil+1);
					elemento[numElementos].nodo[7] = createNode(capa+1, col, fil+1);
					elemento[numElementos].numNodos = 8;
					elemento[numElementos].col = col;
					elemento[numElementos].fil = fil;
					elemento[numElementos].capa = capa;
					numElementos++;
					numDatos += 9;
				}
			}
		}
	}
	return true;
}

//---------------------------------------------------------------------------------
//----------- createNode ----------------------------------------------------------
//---------------------------------------------------------------------------------
int MohidResults::createNode(int capa, int col, int fil)
{
	int capaMask;
	if (capa < numCapas) {
		// Es una capa de celdas
		capaMask = capa;
	} else {
		// Es la capa que está por encima de la superficie libre
		capaMask = capa-1;
	}
	if (nodollbCelda[capa][col][fil] == -1) {
		//El nodo no existe --> crearlo
		double	sumaZ ;
		int		numSumaZ ;
		nodo.resize(nodo.size()+1);
		nodo[numNodos].x = x[col][fil];
		nodo[numNodos].y = y[col][fil];
		nodo[numNodos].capa = capa;
		sumaZ = 0;
		numSumaZ = 0;
		//Calcular la coordenada Z promedio de las celdas vecinas en fias y columnas
		if (col > 0) {
			if (mask[capaMask][col-1][fil] == 1) {
				numSumaZ++;
				sumaZ += z[capa][col-1][fil];
			}
			if (fil > 0) {
				if (mask[capaMask][col-1][fil-1] == 1) {
					numSumaZ++;
					sumaZ += z[capa][col-1][fil-1];
				}
			}
		}
		if (fil > 0) {
			if (mask[capaMask][col][fil-1] == 1) {
				numSumaZ++;
				sumaZ += z[capa][col][fil-1];
			}
		}
		if (mask[capaMask][col][fil] == 1) {
			numSumaZ++;
			sumaZ += z[capa][col][fil];
		}

		nodo[numNodos].z = sumaZ / (double)numSumaZ;

		nodollbCelda[capa][col][fil] = numNodos;
		numNodos++;
		return numNodos-1;
	} else {
		//El nodo existe
		return nodollbCelda[capa][col][fil];
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

	valorNodo.resize(numNodos, 0);
	valorNodo2.resize(numNodos, 0);
	valorNodo3.resize(numNodos, 0);

	cout << "    Escribiendo los resultados en el archivo vtk..." << "\n";

	// open file for output
	ofstream vtk(archVTK);

	//Encabezado
	vtk << "# vtk DataFile Version 2.0" << "\n";
	vtk << archVTK << ", Created by Gmsh" << "\n";
	vtk << "ASCII" << "\n";
	vtk << "DATASET UNSTRUCTURED_GRID" << "\n";

	//Nodos
	vtk << "POINTS " << numNodos << " double" << "\n";
	for (i=0; i<numNodos; i++) {
		if (EscribirComo2D == true) {
			vtk << std::setprecision(10) << nodo[i].x << " " << nodo[i].y << " " << nodo[i].capa * 1.0 << "\n";
		} else {
			vtk << std::setprecision(10) << nodo[i].x << " " << nodo[i].y << " " << nodo[i].z  + offsetVertical << "\n";
		}
	}
	vtk << "\n";

	//Celdas
	vtk << "CELLS " << numElementos << " " << numDatos << "\n";
	for (i=0; i<numElementos; i++) {
		vtk << elemento[i].numNodos << " " << elemento[i].nodo[0]
									<< " " << elemento[i].nodo[1]
									<< " " << elemento[i].nodo[2]
									<< " " << elemento[i].nodo[3]
									<< " " << elemento[i].nodo[4]
									<< " " << elemento[i].nodo[5]
									<< " " << elemento[i].nodo[6]
									<< " " << elemento[i].nodo[7] << "\n";
	}
	vtk << "\n";
	vtk << "CELL_TYPES " << numElementos << "\n";
	for (i=0; i<numElementos; i++) {
		vtk << "12" << "\n"; //VTK_HEXA
	}
	vtk << "\n";

	//Valores en Celdas --------------------------------------
	vtk << "CELL_DATA " << numElementos << "\n";
	//Valores de Velocidad
	vtk << "VECTORS " << "vel" << " float" << "\n";
	for (i=0; i<numElementos; i++) {
		vtk << u[elemento[i].capa][elemento[i].col][elemento[i].fil]
			<< " " << v[elemento[i].capa][elemento[i].col][elemento[i].fil]
			<< " " << w[elemento[i].capa][elemento[i].col][elemento[i].fil] << "\n";
	}
	vtk << "\n";

	//Valores de tirante
	vtk << "SCALARS " << "Tirante" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numElementos; i++) {
		vtk << z[numCapas][elemento[i].col][elemento[i].fil] - z[0][elemento[i].col][elemento[i].fil] << "\n";
	}
	vtk << "\n";

	//Valores de campos adicionales
	for (int nc=0; nc<listaCampos.size(); nc++) {
		vtk << "SCALARS " << listaCampos[nc] << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numElementos; i++) {
			vtk << campos[nc][elemento[i].capa][elemento[i].col][elemento[i].fil] << "\n";
		}
		vtk << "\n";
	}

	//Valores de salinidad
	if (hasSalinity) {
		vtk << "SCALARS " << "salinity" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numElementos; i++) {
			vtk << salinity[elemento[i].capa][elemento[i].col][elemento[i].fil] << "\n";
		}
		vtk << "\n";
	}
	//Valores de densidad
	if (hasDensity) {
		vtk << "SCALARS " << "density" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numElementos; i++) {
			vtk << density[elemento[i].capa][elemento[i].col][elemento[i].fil] << "\n";
		}
		vtk << "\n";
	}
	//Valores de coliformes
	if (hasColiforms) {
		vtk << "SCALARS " << "coliforms" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numElementos; i++) {
			vtk << coliforms[elemento[i].capa][elemento[i].col][elemento[i].fil] << "\n";
		}
		vtk << "\n";
	}
	//Valores de temperatura
	if (hasTemperature) {
		vtk << "SCALARS " << "temperature" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numElementos; i++) {
			vtk << temperature[elemento[i].capa][elemento[i].col][elemento[i].fil] << "\n";
		}
		vtk << "\n";
	}
	//Valores del mapa
	if (cargadoMapa) {
		vtk << "SCALARS " << "Mapa" << " float" << "\n";
		vtk << "LOOKUP_TABLE default" << "\n";
		for (i=0; i<numElementos; i++) {
			vtk << mapa[elemento[i].col][elemento[i].fil] << "\n";
		}
		vtk << "\n";
	}

	// Valores en Nodos --------------------------------------
	vtk << "POINT_DATA " << numNodos << "\n";
	//Valores de batimetria
	calculateNodeValues2D(batimetria, valorNodo);
	vtk << "SCALARS " << "batimetria" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numNodos; i++) {
		vtk << offsetVertical - valorNodo[i] << "\n";
	}

	vtk << "\n";

	//Valores de nivel de agua
	calculateNodeValues2D(nivelSuperficie, valorNodo);
	vtk << "SCALARS " << "NivelSuperficie" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numNodos; i++) {
		vtk << valorNodo[i] + offsetVertical << "\n";
	}

	vtk << "\n";

	//Valores de Velocidad
	calculateNodeValues(u, valorNodo);
	calculateNodeValues(v, valorNodo2);
	calculateNodeValues(w, valorNodo3);
	vtk << "VECTORS " << "vel" << " float" << "\n";
	for (i=0; i<numNodos; i++) {
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

	valorNodo.resize(numNodos, 0);
	valorNodo2.resize(numNodos, 0);
	valorNodo3.resize(numNodos, 0);

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
	vtk << "POINTS " << numNodos << " float" << "\n";

	lastPos = vtk.tellp();
	for (i=0; i<numNodos; i++) {
		dato = nodo[i].x;
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
		dato = nodo[i].y;
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
		dato = nodo[i].z + offsetVertical;
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
		newPos = vtk.tellp();
		if (newPos-lastPos != 12) {
			cout << "    Mal nodo " << i << " - " << lastPos << " " << newPos << "\n";
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
	vtk << "CELLS " << numElementos << " " << numDatos << "\n";
	lastPos = vtk.tellp();
	for (i=0; i<numElementos; i++) {
		datoInt = elemento[i].numNodos;
		vtk.write(reinterpret_cast<const char *>(&datoInt)+3,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt)+2,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt)+1,sizeof(int)/4);
		vtk.write(reinterpret_cast<const char *>(&datoInt),sizeof(int)/4);

		for (j=0; j<8; j++) {
			datoInt = elemento[i].nodo[j];
			vtk.write(reinterpret_cast<const char *>(&datoInt)+3,sizeof(int)/4);
			vtk.write(reinterpret_cast<const char *>(&datoInt)+2,sizeof(int)/4);
			vtk.write(reinterpret_cast<const char *>(&datoInt)+1,sizeof(int)/4);
			vtk.write(reinterpret_cast<const char *>(&datoInt),sizeof(int)/4);
		}
		newPos = vtk.tellp();
		if (newPos-lastPos != 36) {
			cout << "    Mal elemento " << i << " - " << lastPos << " " << newPos << "\n";
			return false;
		}
		lastPos = newPos;
	}
	vtk << "\n";
	vtk << "CELL_TYPES " << numElementos << "\n";
	lastPos = vtk.tellp();
	for (i=0; i<numElementos; i++) {
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
	vtk << "CELL_DATA " << numElementos << "\n";
	//Valores de Velocidad
	vtk << "VECTORS " << "vel" << " float" << "\n";
	for (i=0; i<numElementos; i++) {
		//vtk << u[elemento[i].capa][elemento[i].col][elemento[i].fil]
		//	<< " " << v[elemento[i].capa][elemento[i].col][elemento[i].fil]
		//	<< " " << w[elemento[i].capa][elemento[i].col][elemento[i].fil] << "\n";
		dato = u[elemento[i].capa][elemento[i].col][elemento[i].fil];
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
		dato = v[elemento[i].capa][elemento[i].col][elemento[i].fil];
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
		dato = w[elemento[i].capa][elemento[i].col][elemento[i].fil];
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
	}
	vtk << "\n";

	// Valores en Nodos --------------------------------------
	vtk << "POINT_DATA " << numNodos << "\n";
	//Valores de batimetria
	calculateNodeValues2D(batimetria, valorNodo);
	vtk << "SCALARS " << "batimetria" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numNodos; i++) {
		dato = valorNodo[i] + offsetVertical;
		vtk.write(reinterpret_cast<const char *>(&dato)+3,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+2,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato)+1,sizeof(float)/4);
		vtk.write(reinterpret_cast<const char *>(&dato),sizeof(float)/4);
	}
	vtk << "\n";
	vtk << "\n";

	//Valores de nivel de agua
	calculateNodeValues2D(nivelSuperficie, valorNodo);
	vtk << "SCALARS " << "NivelSuperficie" << " float" << "\n";
	vtk << "LOOKUP_TABLE default" << "\n";
	for (i=0; i<numNodos; i++) {
		dato = valorNodo[i] + offsetVertical;
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
		for (i=0; i<numNodos; i++) {
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
	for (i=0; i<numNodos; i++) {
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

	valorNodo.resize(numNodos, 0);
	valorNodo2.resize(numNodos, 0);
	valorNodo3.resize(numNodos, 0);

	cout << "    Escribiendo los resultados en formato del GIS..." << "\n";
	{
		getDatasetName("supLibre_", indice, nombreArchivo);
		strcat ( nombreArchivo, ".asc");

		cout << "\t\tEscribiendo " << nombreArchivo << "\n";
		std::ofstream surfer(nombreArchivo);

		surfer << "ncols         " << numCol << "\n";
		surfer << "nrows         " << numFil << "\n";
		surfer << "xllcorner     " << x[0][0] << "\n";
		surfer << "yllcorner     " << y[0][0] << "\n";
		surfer << "cellsize      " << x[1][0] - x[0][0]<< "\n";
		surfer << "NODATA_value  -9999.0" << "\n";

		for (int k=numFil-1; k>=0; k--) {
			for (int i=0; i<numCol; i++) {
				if (mask[0][i][k] != 0) {
					surfer << z[numCapas][i][k] + offsetVertical << " ";
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
		surfer << "nrows         " << numFil << "\n";
		surfer << "xllcorner     " << x[0][0] << "\n";
		surfer << "yllcorner     " << y[0][0] << "\n";
		surfer << "cellsize      " << x[1][0] - x[0][0]<< "\n";
		surfer << "NODATA_value  -9999.0" << "\n";

		for (int k=numFil-1; k>=0; k--) {
			for (int i=0; i<numCol; i++) {
				if (mask[0][i][k] != 0) {
					surfer << z[0][i][k] + offsetVertical << " ";
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
		surfer << "nrows         " << numFil << "\n";
		surfer << "xllcorner     " << x[0][0] << "\n";
		surfer << "yllcorner     " << y[0][0] << "\n";
		surfer << "cellsize      " << x[1][0] - x[0][0]<< "\n";
		surfer << "NODATA_value  -9999.0" << "\n";

		for (int k=numFil-1; k>=0; k--) {
			for (int i=0; i<numCol; i++) {
				if (mask[0][i][k] != 0) {
					double um = 0, vm = 0, wm = 0, modU = 0;
					double tirante = 0;
					for (int capa = 0; capa < numCapas; capa++) {
						double dh = z[capa+1][i][k] - z[capa][i][k];
						um += u[capa][i][k] * dh;
						vm += v[capa][i][k] * dh;
						wm += w[capa][i][k] * dh;
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

	for (int nc=0; nc<listaCampos.size(); nc++) {
		getDatasetName(listaCampos[nc].data(), indice, nombreArchivo);
		strcat ( nombreArchivo, ".asc");

		cout << "\t\tEscribiendo " << nombreArchivo << "\n";
		std::ofstream surfer(nombreArchivo);

		surfer << "ncols         " << numCol << "\n";
		surfer << "nrows         " << numFil << "\n";
		surfer << "xllcorner     " << x[0][0] << "\n";
		surfer << "yllcorner     " << y[0][0] << "\n";
		surfer << "cellsize      " << x[1][0] - x[0][0]<< "\n";
		surfer << "NODATA_value  -9999.0" << "\n";

		for (int k=numFil-1; k>=0; k--) {
			for (int i=0; i<numCol; i++) {
				if (mask[0][i][k] != 0) {
					double medio = 0;
					double tirante = 0;
					for (int capa = 0; capa < numCapas; capa++) {
						double dh = z[capa+1][i][k] - z[capa][i][k];
						medio += campos[nc][capa][i][k] * dh;
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

	for (i=0; i<numNodos; i++) {
		valNodo[i]=0;
	}

	for (i=0; i<numElementos; i++) {
		for (j=0; j<elemento[i].numNodos; j++) {
			valNodo[elemento[i].nodo[j]] += valElemento[elemento[i].capa][elemento[i].col][elemento[i].fil];
			numValNodo[elemento[i].nodo[j]]++;
		}
	}
	for (i=0; i<numNodos; i++) {
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

	for (i=0; i<numNodos; i++) {
		valNodo[i]=0;
	}

	for (i=0; i<numElementos; i++) {
		for (j=0; j<elemento[i].numNodos; j++) {
			valNodo[elemento[i].nodo[j]] += valElemento[elemento[i].col][elemento[i].fil];
			numValNodo[elemento[i].nodo[j]]++;
		}
	}
	for (i=0; i<numNodos; i++) {
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
