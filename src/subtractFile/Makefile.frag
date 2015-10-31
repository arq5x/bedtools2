BUILT_OBJECTS += obj/subtractHelp.o obj/subtractFile.o

obj/subtractHelp.o: src/subtractFile/subtractHelp.cpp obj/subtractHelp.d
	$(CXX_COMPILE)

obj/subtractFile.o: src/subtractFile/subtractFile.cpp obj/subtractFile.d
	$(CXX_COMPILE)
