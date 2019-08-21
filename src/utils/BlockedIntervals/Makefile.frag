BUILT_OBJECTS += obj/BlockedIntervals.o

obj/BlockedIntervals.o: src/utils/BlockedIntervals/BlockedIntervals.cpp obj/BlockedIntervals.d
	$(CXX_COMPILE)
