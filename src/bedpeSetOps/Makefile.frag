BUILT_OBJECTS += obj/bedpeSetOps.o

obj/bedpeSetOps.o: src/bedpeSetOps/bedpeSetOps.cpp obj/bedpeSetOps.d
	$(CXX_COMPILE)
