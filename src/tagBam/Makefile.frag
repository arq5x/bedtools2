BUILT_OBJECTS += obj/tagBamMain.o obj/tagBam.o

obj/tagBamMain.o: src/tagBam/tagBamMain.cpp obj/tagBamMain.d
	$(CXX_COMPILE)

obj/tagBam.o: src/tagBam/tagBam.cpp obj/tagBam.d
	$(CXX_COMPILE)
