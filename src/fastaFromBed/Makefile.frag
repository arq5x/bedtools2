BUILT_OBJECTS += obj/fastaFromBedMain.o obj/fastaFromBed.o

obj/fastaFromBedMain.o: src/fastaFromBed/fastaFromBedMain.cpp obj/fastaFromBedMain.d
	$(CXX_COMPILE)

obj/fastaFromBed.o: src/fastaFromBed/fastaFromBed.cpp obj/fastaFromBed.d
	$(CXX_COMPILE)
