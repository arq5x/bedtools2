BUILT_OBJECTS += obj/linksMain.o obj/linksBed.o

obj/linksMain.o: src/linksBed/linksMain.cpp obj/linksMain.d
	$(CXX_COMPILE)

obj/linksBed.o: src/linksBed/linksBed.cpp obj/linksBed.d
	$(CXX_COMPILE)
