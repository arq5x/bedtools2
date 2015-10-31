BUILT_OBJECTS += obj/closestHelp.o obj/closestFile.o

obj/closestHelp.o: src/closestFile/closestHelp.cpp obj/closestHelp.d
	$(CXX_COMPILE)

obj/closestFile.o: src/closestFile/closestFile.cpp obj/closestFile.d
	$(CXX_COMPILE)
