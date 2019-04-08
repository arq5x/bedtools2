BUILT_OBJECTS += obj/NewChromsweep.o obj/CloseSweep.o

obj/NewChromsweep.o: src/utils/NewChromsweep/NewChromsweep.cpp obj/NewChromsweep.d
	$(CXX_COMPILE)

obj/CloseSweep.o: src/utils/NewChromsweep/CloseSweep.cpp obj/CloseSweep.d
	$(CXX_COMPILE)
