BUILT_OBJECTS += obj/slopBedMain.o obj/slopBed.o

obj/slopBedMain.o: src/slopBed/slopBedMain.cpp obj/slopBedMain.d
	$(CXX_COMPILE)

obj/slopBed.o: src/slopBed/slopBed.cpp obj/slopBed.d
	$(CXX_COMPILE)
