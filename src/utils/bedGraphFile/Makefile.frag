BUILT_OBJECTS += obj/bedGraphFile.o

obj/bedGraphFile.o: src/utils/bedGraphFile/bedGraphFile.cpp obj/bedGraphFile.d
	$(CXX_COMPILE)
