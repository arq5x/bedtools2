BUILT_OBJECTS += obj/sampleHelp.o obj/sampleFile.o

obj/sampleHelp.o: src/sampleFile/sampleHelp.cpp obj/sampleHelp.d
	$(CXX_COMPILE)

obj/sampleFile.o: src/sampleFile/sampleFile.cpp obj/sampleFile.d
	$(CXX_COMPILE)
