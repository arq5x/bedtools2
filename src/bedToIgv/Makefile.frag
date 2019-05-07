BUILT_OBJECTS += obj/bedToIgv.o

obj/bedToIgv.o: src/bedToIgv/bedToIgv.cpp obj/bedToIgv.d
	$(CXX_COMPILE)
