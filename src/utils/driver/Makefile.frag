BUILT_OBJECTS += obj/BedtoolsDriver.o

obj/BedtoolsDriver.o: src/utils/driver/BedtoolsDriver.cpp obj/BedtoolsDriver.d
	$(CXX_COMPILE)
