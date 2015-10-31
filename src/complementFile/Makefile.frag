BUILT_OBJECTS += obj/complementHelp.o obj/complementFile.o

obj/complementHelp.o: src/complementFile/complementHelp.cpp obj/complementHelp.d
	$(CXX_COMPILE)

obj/complementFile.o: src/complementFile/complementFile.cpp obj/complementFile.d
	$(CXX_COMPILE)
