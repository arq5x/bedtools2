BUILT_OBJECTS += obj/windowMain.o obj/windowBed.o

obj/windowMain.o: src/windowBed/windowMain.cpp obj/windowMain.d
	$(CXX_COMPILE)

obj/windowBed.o: src/windowBed/windowBed.cpp obj/windowBed.d
	$(CXX_COMPILE)
