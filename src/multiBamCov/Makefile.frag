BUILT_OBJECTS += obj/multiBamCovMain.o obj/multiBamCov.o

obj/multiBamCovMain.o: src/multiBamCov/multiBamCovMain.cpp obj/multiBamCovMain.d
	$(CXX_COMPILE)

obj/multiBamCov.o: src/multiBamCov/multiBamCov.cpp obj/multiBamCov.d
	$(CXX_COMPILE)
