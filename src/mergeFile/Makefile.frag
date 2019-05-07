BUILT_OBJECTS += obj/mergeHelp.o obj/mergeFile.o

obj/mergeHelp.o: src/mergeFile/mergeHelp.cpp obj/mergeHelp.d
	$(CXX_COMPILE)

obj/mergeFile.o: src/mergeFile/mergeFile.cpp obj/mergeFile.d
	$(CXX_COMPILE)
