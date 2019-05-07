BUILT_OBJECTS += obj/splitBedMain.o obj/splitBed.o

obj/splitBedMain.o: src/split/splitBedMain.cpp obj/splitBedMain.d
	$(CXX_COMPILE)

obj/splitBed.o: src/split/splitBed.cpp obj/splitBed.d
	$(CXX_COMPILE)
