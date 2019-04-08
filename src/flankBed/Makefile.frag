BUILT_OBJECTS += obj/flankBedMain.o obj/flankBed.o

obj/flankBedMain.o: src/flankBed/flankBedMain.cpp obj/flankBedMain.d
	$(CXX_COMPILE)

obj/flankBed.o: src/flankBed/flankBed.cpp obj/flankBed.d
	$(CXX_COMPILE)
