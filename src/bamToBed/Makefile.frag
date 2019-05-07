BUILT_OBJECTS += obj/bamToBed.o

obj/bamToBed.o: src/bamToBed/bamToBed.cpp obj/bamToBed.d
	$(CXX_COMPILE)
