BUILT_OBJECTS += obj/genomeCoverageMain.o obj/genomeCoverageBed.o

obj/genomeCoverageMain.o: src/genomeCoverageBed/genomeCoverageMain.cpp obj/genomeCoverageMain.d
	$(CXX_COMPILE)

obj/genomeCoverageBed.o: src/genomeCoverageBed/genomeCoverageBed.cpp obj/genomeCoverageBed.d
	$(CXX_COMPILE)
