# Define the symbols we might want to change:
CXX	:= gcc
CXXFLAGS := -g

OBJECTS	:= MohidResult.o  mohidToVTK.o

#-L/usr
LIBS	:= -lstdc++ -lhdf5 -lhdf5_hl -lz
LIBSDIR	:= -L../../../hdf5-1.8.4-patch1/hdf5/lib -L../../../zlib-1.2.5
INCLUDES := -I../../../hdf5-1.8.4-patch1/hdf5/include

mohidToVTK: $(OBJECTS)
	$(CXX) $(OBJECTS) -o mohidToVTK $(LIBS) $(LIBSDIR)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
