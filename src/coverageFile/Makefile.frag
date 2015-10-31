BUILT_OBJECTS += obj/coverageHelp.o obj/coverageFile.o

obj/coverageHelp.o: src/coverageFile/coverageHelp.cpp obj/coverageHelp.d
	$(CXX_COMPILE)

obj/coverageFile.o: src/coverageFile/coverageFile.cpp obj/coverageFile.d
	$(CXX_COMPILE)
