BUILT_OBJECTS += obj/nucBedMain.o obj/nucBed.o

obj/nucBedMain.o: src/nucBed/nucBedMain.cpp obj/nucBedMain.d
	$(CXX_COMPILE)

obj/nucBed.o: src/nucBed/nucBed.cpp obj/nucBed.d
	$(CXX_COMPILE)
