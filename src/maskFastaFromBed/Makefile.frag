BUILT_OBJECTS += obj/maskFastaFromBedMain.o obj/maskFastaFromBed.o

obj/maskFastaFromBedMain.o: src/maskFastaFromBed/maskFastaFromBedMain.cpp obj/maskFastaFromBedMain.d
	$(CXX_COMPILE)

obj/maskFastaFromBed.o: src/maskFastaFromBed/maskFastaFromBed.cpp obj/maskFastaFromBed.d
	$(CXX_COMPILE)
