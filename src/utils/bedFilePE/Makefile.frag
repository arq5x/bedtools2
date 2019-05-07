BUILT_OBJECTS += obj/bedFilePE.o

obj/bedFilePE.o: src/utils/bedFilePE/bedFilePE.cpp obj/bedFilePE.d
	$(CXX_COMPILE)
