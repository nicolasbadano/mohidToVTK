**mohidToVTK** is a tool to convert MOHID results into the VTK and Raster formats. It is written in C++ and requires for compiling the HDF5, ZLIB and SLIB libraries.

mohidToVTK is entirely **open-source**, and everyone is free to propose changes or maintain their own, customized version as long as they make their changes open to the public in accordance with the GNU General Public License (for more information check the license file attached to this project).

### Current Features
*   Converts MOHID results to 2D or 3D VTK files.
*   Converts MOHID results to ASC Raster format.
*   Can retrieve any result field, including Hydrodynamic, WaterProperties and Lagrangian fields.

### Running
Just run the compiled executable file from a console terminal in the folder containing the HDF5 result files. To convert default hydrodynamic results, just run:

`mohidToVTK Hydrodynamic_1.hdf5`

To convert specific fields, specify a file containing the required field paths using the `-c` flag, and a list of source files:

`mohidToVTK -c fields.txt Hydrodynamic_1.hdf5 Lagrangian_3.hdf5`

An example fields file is provided on the `testdata` folder. The field path for any given field can be obtained inspecting the results file with a HDF5 inspector program. For water quality parameters, they typically take the following format:

```
/Results/temperature/temperature
/Results/fecal coliforms/fecal coliforms
/Results/salinity/salinity
/Results/density/density
```

Meanwhile, lagrangian fields typically look like:

```
/Results/Group_1/Data_3D/fecal coliforms/fecal coliforms
```

#### Compiling
A project file for Visual Studio 2010 is provided for compilation under Windows. To compile under linux, a simple `makefile` is available on the repository. In both cases paths need to be provided both for the  include files and the binaries of the required libraries: hdf5, zlib and slib.

### Contributing
If you want to help put with the ongoing development, you can do so by looking for possible bugs or by contributing new features. To become one of the developers, simply fork this repository and submit your pull requests for review.

To report a bug, propose a feature, or suggest a change to the existing code please, use our [Issue Tracker](https://github.com/esteldunedain/mohidtovtk/issues).

### Author
Written by Nicolás Diego Badano < nicolas.d.badano at gmail.com >, National Institute for Water, Argentina
