BUILT_OBJECTS += obj/sortMain.o obj/sortBed.o

obj/sortMain.o: src/sortBed/sortMain.cpp obj/sortMain.d
	$(CXX_COMPILE)

obj/sortBed.o: src/sortBed/sortBed.cpp obj/sortBed.d
	$(CXX_COMPILE)
