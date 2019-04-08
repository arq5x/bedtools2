BUILT_OBJECTS += obj/getOverlap.o

obj/getOverlap.o: src/getOverlap/getOverlap.cpp obj/getOverlap.d
	$(CXX_COMPILE)
