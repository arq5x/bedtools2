BUILT_OBJECTS += obj/bedpeToBam.o

obj/bedpeToBam.o: src/bedpeToBam/bedpeToBam.cpp obj/bedpeToBam.d
	$(CXX_COMPILE)
