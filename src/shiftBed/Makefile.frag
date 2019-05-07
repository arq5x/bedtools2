BUILT_OBJECTS += obj/shiftBedMain.o obj/shiftBed.o

obj/shiftBedMain.o: src/shiftBed/shiftBedMain.cpp obj/shiftBedMain.d
	$(CXX_COMPILE)

obj/shiftBed.o: src/shiftBed/shiftBed.cpp obj/shiftBed.d
	$(CXX_COMPILE)
