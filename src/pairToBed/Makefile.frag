BUILT_OBJECTS += obj/pairToBedMain.o obj/pairToBed.o

obj/pairToBedMain.o: src/pairToBed/pairToBedMain.cpp obj/pairToBedMain.d
	$(CXX_COMPILE)

obj/pairToBed.o: src/pairToBed/pairToBed.cpp obj/pairToBed.d
	$(CXX_COMPILE)
