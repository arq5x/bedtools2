BUILT_OBJECTS += obj/shuffleBedMain.o obj/shuffleBed.o

obj/shuffleBedMain.o: src/shuffleBed/shuffleBedMain.cpp obj/shuffleBedMain.d
	$(CXX_COMPILE)

obj/shuffleBed.o: src/shuffleBed/shuffleBed.cpp obj/shuffleBed.d
	$(CXX_COMPILE)
