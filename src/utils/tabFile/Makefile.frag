BUILT_OBJECTS += obj/tabFile.o

obj/tabFile.o: src/utils/tabFile/tabFile.cpp obj/tabFile.d
	$(CXX_COMPILE)
