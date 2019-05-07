BUILT_OBJECTS += obj/sequenceUtils.o

obj/sequenceUtils.o: src/utils/sequenceUtilities/sequenceUtils.cpp obj/sequenceUtils.d
	$(CXX_COMPILE)
