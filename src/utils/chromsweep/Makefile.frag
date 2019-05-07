BUILT_OBJECTS += obj/chromsweep.o

obj/chromsweep.o: src/utils/chromsweep/chromsweep.cpp obj/chromsweep.d
	$(CXX_COMPILE)
