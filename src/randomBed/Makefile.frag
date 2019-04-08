BUILT_OBJECTS += obj/randomBedMain.o obj/randomBed.o

obj/randomBedMain.o: src/randomBed/randomBedMain.cpp obj/randomBedMain.d
	$(CXX_COMPILE)

obj/randomBed.o: src/randomBed/randomBed.cpp obj/randomBed.d
	$(CXX_COMPILE)
