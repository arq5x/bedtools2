BUILT_OBJECTS += obj/unionBedGraphsMain.o obj/unionBedGraphs.o

obj/unionBedGraphsMain.o: src/unionBedGraphs/unionBedGraphsMain.cpp obj/unionBedGraphsMain.d
	$(CXX_COMPILE)

obj/unionBedGraphs.o: src/unionBedGraphs/unionBedGraphs.cpp obj/unionBedGraphs.d
	$(CXX_COMPILE)
