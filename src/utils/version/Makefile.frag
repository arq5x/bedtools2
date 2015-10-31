BUILT_OBJECTS += obj/version.o

obj/version.o: src/utils/version/version.cpp obj/version.d
	$(CXX_COMPILE)
