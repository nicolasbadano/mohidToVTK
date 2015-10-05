// -----------------------------------------------------------------------
// -------------- mohidToVTK - v2 ----------------------------------------
// -----------------------------------------------------------------------
// -- Converts Mohid Results into VTK ascii format -----------------------
// ----------- by Nicolas Badano --- 01/2013 -----------------------------
// -----------------------------------------------------------------------
// ---------------- A.M.D.G. ---------------------------------------------
// -----------------------------------------------------------------------

#include "ClassMohidResult.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>

int main(int argc, char* argv[])
{
	ClassMohidResult			*hbRes;
	bool						existe;
	int							indice;
	int							numArchivos = 0;
	char						nombreArchivoHDF5[500];
	char						nombreArchivoHDF5_2[500];
	char						nombreArchivoMapa[500];	
	char						nombreArchivoCampos[300];
	bool						haySegundoArchivo;
	char						nombreArchivo[300];
	char						pasos[100];
	int							count;
	int							i, j;
	int							pasoInicial, pasoFinal, cadaCuantos;
	int 						inicio, fin;
	bool						terminar = false;
	bool						cargarMapa = false;
	bool						archivoCampos = false;
	double						offset = 0;

	cout << "\n";
	cout << "----------------------------- mohid3DToVTK3D -----------------------------\n";
	cout << "\n";

	pasoInicial = 1; pasoFinal = 1000000; cadaCuantos = 1;

    // Evaluar parámetros recibidos
	for (i=1; i < argc; i++)
	{
		if (argv[i][0] == '-' && argv[i][1] == 'm')
		{
			strcpy(nombreArchivoMapa, argv[i+1]);
			cargarMapa = true;
			i++;
			continue;
		}
		if (argv[i][0] == '-' && argv[i][1] == 'c')
		{
			strcpy(nombreArchivoCampos, argv[i+1]);
			archivoCampos = true;
			i++;
			continue;
		}		
		if (argv[i][0] == '-' && argv[i][1] == 'o')
		{
			offset = strtod(argv[i+1], 0);
			i++;
			continue;
		}
		
		if (argv[i][0] == '-' && argv[i][1] == 't')
		{
			// Se está especificando paso inicial:pasofinal:cadacuantos
			inicio = 2;
			for (j = inicio; terminar == false; j++)
			{
				if (argv[i][j] == ':' || argv[i][j] == '\0')
				{
					if (argv[i][j] == '\0') terminar = true;
					fin = j;
					if (fin > inicio)
					{
						argv[i][j] = '\0';						
						// Encontramos el primer valor, el paso inicial
						pasoInicial = atoi(argv[i]+inicio);		
						pasoFinal = pasoInicial;
					}
					inicio = fin + 1;	
					break;
				}
			}
			for (j = inicio; terminar == false; j++)
			{
				if (argv[i][j] == ':' || argv[i][j] == '\0')
				{
					if (argv[i][j] == '\0') terminar = true;
					fin = j;
					if (fin > inicio)
					{
						// Encontramos el segundo valor, el paso final
						argv[i][j] = '\0';
						pasoFinal = atoi(argv[i]+inicio);		
					}
					inicio = fin + 1;	
					break;
				}
			}
			for (j = inicio; terminar == false; j++)
			{
				if (argv[i][j] == '\0')
				{
					terminar = true;
					fin = j;
					if (fin > inicio)
					{
						// Encontramos el tercer valor, cada cuantos
						argv[i][j] = '\0';
						cadaCuantos = atoi(argv[i]+inicio);		
					}
					inicio = fin + 1;	
					break;
				}
			}
		}
		else
		{
			// Asumimos que es un nombre de archivo
			if (numArchivos == 0)
			{
				strcpy(nombreArchivoHDF5, argv[i]);
				numArchivos++;
			}
			else if (numArchivos == 1)
			{
				strcpy(nombreArchivoHDF5_2, argv[i]);
				numArchivos++;
			} else ; //Ignorar archivos subsiguientes
		}
	}
	if (numArchivos == 0)
	{	// No pasaron nombre de archivo
		strcpy(nombreArchivoHDF5, "Hydrodynamic_1.hdf5");
		numArchivos = 1;
	}

	cout << "Paso Inicial a Exportar: " << pasoInicial << "\n";
	cout << "Paso Final a Exportar: " << pasoFinal << "\n";
	cout << "Cada Cuantos pasos Exportar: " << cadaCuantos << "\n";

	existe = true;
	indice = pasoInicial;
	while (existe && indice <= pasoFinal)
	{
		hbRes = new ClassMohidResult(offset);

		existe = hbRes->CargarResultado(nombreArchivoHDF5, indice);
		if (existe)
		{
			if (numArchivos == 2)
			{			
				if (archivoCampos)
					existe = hbRes->cargarCamposArchivo( nombreArchivoHDF5_2, nombreArchivoCampos, indice);
				else
					existe = hbRes->AgregarResultado(nombreArchivoHDF5_2, indice);

				if (existe == false) break;
			}

			if (cargarMapa)
				hbRes->cargarMapa( nombreArchivoMapa );

			//Convertir los resultados a VTK
			hbRes->ConvertirResultadoAVTK();

			//Escribir el archivo VTK
			hbRes->DevolverNombreDataset("Hydro_", indice, nombreArchivo);
			strcat(nombreArchivo, ".vtk");
			hbRes->EscribirResultadoVTK(nombreArchivo,false);
			
			hbRes->DevolverNombreDataset("Hydro2D_", indice, nombreArchivo);
			strcat(nombreArchivo, ".vtk");
			hbRes->EscribirResultadoVTK(nombreArchivo,true);
			
			hbRes->EscribirResultadoGIS(indice);
		}
	
		delete hbRes;
		indice+=cadaCuantos;
	}
	return 0;
}

