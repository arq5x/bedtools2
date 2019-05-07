BUILT_OBJECTS += obj/bedFile.o

obj/bedFile.o: src/utils/bedFile/bedFile.cpp obj/bedFile.d
	$(CXX_COMPILE)
