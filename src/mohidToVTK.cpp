// mohidToVTK, the utility to convert MOHID results to VTK and Raster formats-
// Copyright (C) 2010-2015  Nicolás Diego Badano
// GNU GPLv2
// (c) Nicolás Diego Badano,
// A.M.D.G.

#include "MohidResults.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>

int main(int argc, char* argv[])
{
    cout << "\n";
    cout << "----------------------------- mohidToVTK -----------------------------\n";
    cout << "\n";

    int numFiles = 0;
    char hdf5FileName1[500];
    char hdf5FileName2[500];
    char fieldsFileName[300];
    char mapFileNane[500];

    bool finish = false;
    bool loadMap = false;
    bool fieldsFileExists = false;
    double offset = 0;

    int stepInitial = 1;
    int stepFinal = 1000000;
    int stepInterval = 1;

    // Evaluate command line parameters
    for (int i=1; i < argc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'm') {
            strcpy(mapFileNane, argv[i+1]);
            loadMap = true;
            i++;
            continue;
        }
        if (argv[i][0] == '-' && argv[i][1] == 'c') {
            strcpy(fieldsFileName, argv[i+1]);
            fieldsFileExists = true;
            i++;
            continue;
        }
        if (argv[i][0] == '-' && argv[i][1] == 'o') {
            offset = strtod(argv[i+1], 0);
            i++;
            continue;
        }

        if (argv[i][0] == '-' && argv[i][1] == 't') {
            // Found initialStep:stepFinal:stepInterval
            int start = 2;
            for (int j = start; finish == false; j++) {
                if (argv[i][j] == ':' || argv[i][j] == '\0') {
                    if (argv[i][j] == '\0') finish = true;
                    int end = j;
                    if (end > start) {
                        argv[i][j] = '\0';
                        // Found the first value, initial step
                        stepInitial = atoi(argv[i]+start);
                        stepFinal = stepInitial;
                    }
                    start = end + 1;
                    break;
                }
            }
            for (int j = start; finish == false; j++) {
                if (argv[i][j] == ':' || argv[i][j] == '\0') {
                    if (argv[i][j] == '\0') finish = true;
                    int end = j;
                    if (end > start) {
                        // Found the second value, final step
                        argv[i][j] = '\0';
                        stepFinal = atoi(argv[i]+start);
                    }
                    start = end + 1;
                    break;
                }
            }
            for (int j = start; finish == false; j++) {
                if (argv[i][j] == '\0') {
                    finish = true;
                    int end = j;
                    if (end > start) {
                        // Found the third value, step interval
                        argv[i][j] = '\0';
                        stepInterval = atoi(argv[i]+start);
                    }
                    start = end + 1;
                    break;
                }
            }
        } else {
            // Asume it's a filename
            if (numFiles == 0) {
                strcpy(hdf5FileName1, argv[i]);
                numFiles++;
            } else if (numFiles == 1) {
                strcpy(hdf5FileName2, argv[i]);
                numFiles++;
            }
            // Ignore subsequent files
        }
    }
    // If no file name was supplied, addopt a default
    if (numFiles == 0) {
        strcpy(hdf5FileName1, "Hydrodynamic_1.hdf5");
        numFiles = 1;
    }

    cout << "Exporting steps " << stepInitial << " through " << stepFinal << " every " << stepInterval << "\n";

    bool exists = true;
    int step = stepInitial;
    while (exists && step <= stepFinal) {
        MohidResults mohidResults(offset);

        exists = mohidResults.loadResult(hdf5FileName1, step);
        if (exists) {
            if (numFiles == 2) {
                if (fieldsFileExists) {
                    exists = mohidResults.loadFieldsFromFieldFile( hdf5FileName2, fieldsFileName, step);
                } else {
                    exists = mohidResults.addResult(hdf5FileName2, step);
                }

                if (exists == false) break;
            }

            if (loadMap) {
                mohidResults.loadMap( mapFileNane );
            }

            //Convertir los resultados a VTK
            mohidResults.convertResultsToVTK();

            //Escribir el archivo VTK
            char fileName[300];
            mohidResults.getDatasetName("Hydro_", step, fileName);
            strcat(fileName, ".vtk");
            mohidResults.writeResultsVTK(fileName, false);

            mohidResults.getDatasetName("Hydro2D_", step, fileName);
            strcat(fileName, ".vtk");
            mohidResults.writeResultsVTK(fileName, true);

            mohidResults.writeResultsASC(step);
        }

        step += stepInterval;
    }
    return 0;
}

