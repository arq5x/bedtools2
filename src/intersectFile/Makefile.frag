BUILT_OBJECTS += obj/intersectHelp.o obj/intersectFile.o

obj/intersectHelp.o: src/intersectFile/intersectHelp.cpp obj/intersectHelp.d
	$(CXX_COMPILE)

obj/intersectFile.o: src/intersectFile/intersectFile.cpp obj/intersectFile.d
	$(CXX_COMPILE)
