BUILT_OBJECTS += obj/pairToPairMain.o obj/pairToPair.o

obj/pairToPairMain.o: src/pairToPair/pairToPairMain.cpp obj/pairToPairMain.d
	$(CXX_COMPILE)

obj/pairToPair.o: src/pairToPair/pairToPair.cpp obj/pairToPair.d
	$(CXX_COMPILE)
