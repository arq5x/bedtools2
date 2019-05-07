BUILT_OBJECTS += obj/clusterMain.o obj/clusterBed.o

obj/clusterMain.o: src/clusterBed/clusterMain.cpp obj/clusterMain.d
	$(CXX_COMPILE)

obj/clusterBed.o: src/clusterBed/clusterBed.cpp obj/clusterBed.d
	$(CXX_COMPILE)
