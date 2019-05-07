BUILT_OBJECTS += obj/bedToBam.o

obj/bedToBam.o: src/bedToBam/bedToBam.cpp obj/bedToBam.d
	$(CXX_COMPILE)
