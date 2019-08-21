BUILT_OBJECTS += obj/bed12ToBed6.o

obj/bed12ToBed6.o: src/bed12ToBed6/bed12ToBed6.cpp obj/bed12ToBed6.d
	$(CXX_COMPILE)
