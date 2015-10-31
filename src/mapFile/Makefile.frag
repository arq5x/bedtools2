BUILT_OBJECTS += obj/mapHelp.o obj/mapFile.o

obj/mapHelp.o: src/mapFile/mapHelp.cpp obj/mapHelp.d
	$(CXX_COMPILE)

obj/mapFile.o: src/mapFile/mapFile.cpp obj/mapFile.d
	$(CXX_COMPILE)
