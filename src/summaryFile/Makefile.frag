BUILT_OBJECTS += obj/summaryHelp.o obj/summaryFile.o

obj/summaryHelp.o: src/summaryFile/summaryHelp.cpp obj/summaryHelp.d
	$(CXX_COMPILE)

obj/summaryFile.o: src/summaryFile/summaryFile.cpp obj/summaryFile.d
	$(CXX_COMPILE)
