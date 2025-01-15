BUILT_OBJECTS += obj/revcompBedMain.o obj/revcompBed.o

obj/revcompBedMain.o: src/revcompBed/revcompBedMain.cpp obj/revcompBedMain.d
	$(CXX_COMPILE)

obj/revcompBed.o: src/revcompBed/revcompBed.cpp obj/revcompBed.d
	$(CXX_COMPILE)
