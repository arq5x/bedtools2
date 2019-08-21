BUILT_OBJECTS += obj/multiIntersectBedMain.o obj/multiIntersectBed.o

obj/multiIntersectBedMain.o: src/multiIntersectBed/multiIntersectBedMain.cpp obj/multiIntersectBedMain.d
	$(CXX_COMPILE)

obj/multiIntersectBed.o: src/multiIntersectBed/multiIntersectBed.cpp obj/multiIntersectBed.d
	$(CXX_COMPILE)
