#pragma once

#include <vector>
#include <string>
#include "hdf5.h"
#include "hdf5_hl.h"

using namespace std;

class Node
{
public:
	double		x;
	double		y;
	double		z;
	int			capa;

	double		u;
	double		v;
	double		w;

	Node(void) {};
	~Node(void) {};
};

class Element
{
public:
	int			nodo[8];
	int			numNodos;

	int			col;
	int			fil;
	int			capa;

	Element(void) {};
	~Element(void) {};
};

class MohidResults
{
public:

	int			numCol;
	int			numFil;
	int			numCapas;

	bool		cargadaMalla;
	bool		cargadoResultado;
	bool		cargadoMapa;

	vector<vector<vector<int> > >			mask;
	vector<vector<vector<double> > >		u;
	vector<vector<vector<double> > >		v;
	vector<vector<vector<double> > >		w;
	vector<vector<double> >					batimetria;
	vector<vector<double> >					nivelSuperficie;
	vector<vector<double> >					mapa;

	vector<vector<vector<double> > >		salinity;
	bool									hasSalinity;

	vector<vector<vector<double> > >		density;
	bool									hasDensity;

	vector<vector<vector<double> > >		coliforms;
	bool									hasColiforms;

	vector<vector<vector<double> > >		temperature;
	bool									hasTemperature;

	vector<vector<vector<vector<double> > > > campos;
	vector<std::string >					listaCampos;


	vector<vector<double> >					x;				// Coordenada x del nodo ll de cada fila y columna
	vector<vector<double> >					y;				// Coordenada y del nodo ll de cada fila y columna
	vector<vector<vector<double> > >		z;				// Coordenada z del piso de cada celda (fila y columna y capas)

	vector<vector<vector<int> > >			nodollbCelda;	//Numero de nodo en la parte inferior izquierda anterior de la celda [capa][col][fil]

	vector<Node>							nodo;
	vector<Element>							elemento;
	int										numNodos;
	int										numElementos;

	int										numDatos;
	float									offsetVertical;
	MohidResults(double offset);
	~MohidResults(void);

	bool CargarResultado(char* archHDF5, int indice);
	bool AgregarResultado(char* archHDF5, int indice);
	bool cargarCamposArchivo(char* archH5, char* archCampos, int indice);

	bool CargarCampo3D(hid_t file_id, char* nombreDataset, vector<vector<vector<double> > > &vectorRes);
	bool CargarCampo2D(hid_t file_id, char* nombreDataset, vector<vector<double> > &vectorRes);
	bool cargarMapa(char* archMapa);
	bool ConvertirResultadoAVTK(void);
	void CalcularValoresEnNodos(vector<vector<vector<double> > > &valElemento, vector<double> &valNodo);
	void CalcularValoresEnNodos2D(vector<vector<double> > &valElemento, vector<double> &valNodo);

	void DevolverNombreDataset(const char *nombreBase, int indice, char *nombreFinal);

	int  CrearNodo(int capa, int col, int fil);

	bool EscribirResultadoVTK(char* archVTK, bool EscribirComo2D);
	bool EscribirResultadoVTKBinario(char* archVTK, bool EscribirComo2D);
	bool EscribirResultadoGIS(int indice);
};
