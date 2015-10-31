BUILT_OBJECTS += obj/VectorOps.o

obj/VectorOps.o: src/utils/VectorOps/VectorOps.cpp obj/VectorOps.d
	$(CXX_COMPILE)
