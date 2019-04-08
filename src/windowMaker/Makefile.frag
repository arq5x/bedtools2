BUILT_OBJECTS += obj/windowMakerMain.o obj/windowMaker.o

obj/windowMakerMain.o: src/windowMaker/windowMakerMain.cpp obj/windowMakerMain.d
	$(CXX_COMPILE)

obj/windowMaker.o: src/windowMaker/windowMaker.cpp obj/windowMaker.d
	$(CXX_COMPILE)
