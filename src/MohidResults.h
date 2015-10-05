#pragma once

#include <vector>
#include <string>
#include "hdf5.h"
#include "hdf5_hl.h"

using namespace std;

class Node
{
public:
    double      x;
    double      y;
    double      z;
    int         layer;

    double      u;
    double      v;
    double      w;

    Node(void) {};
    ~Node(void) {};
};

class Cell
{
public:
    int         node[8];
    int         numNodes;

    int         col;
    int         row;
    int         layer;

    Cell(void) {};
    ~Cell(void) {};
};

class MohidResults
{
public:

    int                                     numCol;
    int                                     numRow;
    int                                     numLayers;

    bool                                    mapIsLoaded;

    vector<vector<vector<int> > >           mask;
    vector<vector<vector<double> > >        u;
    vector<vector<vector<double> > >        v;
    vector<vector<vector<double> > >        w;
    vector<vector<double> >                 bathymetry;
    vector<vector<double> >                 surfaceElevation;
    vector<vector<double> >                 map;

    vector<vector<vector<double> > >        salinity;
    bool                                    hasSalinity;

    vector<vector<vector<double> > >        density;
    bool                                    hasDensity;

    vector<vector<vector<double> > >        coliforms;
    bool                                    hasColiforms;

    vector<vector<vector<double> > >        temperature;
    bool                                    hasTemperature;

    vector<vector<vector<vector<double> > > > fields;
    vector<std::string >                    fieldNames;

    // Coordinate x of the ll node of every cell (row and column)
    vector<vector<double> >                 x;
    // Coordinate y of the ll node of every cell (row and column)
    vector<vector<double> >                 y;
    // Coordinate z of the floor of every cell (row, column and layer)
    vector<vector<vector<double> > >        z;
    // Index of the lower, left, bottom node of the cell [layer][col][row]
    vector<vector<vector<int> > >           cellLLBNodeIndex;

    vector<Node>                            node;
    vector<Cell>                            cell;
    int                                     numNodes;
    int                                     numCells;

    int                                     numValues;
    float                                   verticalOffset;
    MohidResults(double offset);
    ~MohidResults(void);

    bool loadResult(char* archHDF5, int indice);
    bool addResult(char* archHDF5, int indice);
    bool loadFieldsFromFieldFile(char* archH5, char* archCampos, int indice);

    bool loadField3D(hid_t file_id, char* nombreDataset, vector<vector<vector<double> > > &vectorRes);
    bool loadField2D(hid_t file_id, char* nombreDataset, vector<vector<double> > &vectorRes);
    bool loadMap(char* archMapa);
    bool convertResultsToVTK(void);
    void calculateNodeValues(vector<vector<vector<double> > > &valElemento, vector<double> &valNodo);
    void calculateNodeValues2D(vector<vector<double> > &valElemento, vector<double> &valNodo);

    void getDatasetName(const char *nombreBase, int indice, char *nombreFinal);

    int  createNode(int layer, int col, int row);

    bool writeResultsVTK(char* archVTK, bool EscribirComo2D);
    bool writeResultsVTKBinary(char* archVTK, bool EscribirComo2D);
    bool writeResultsASC(int indice);
};
