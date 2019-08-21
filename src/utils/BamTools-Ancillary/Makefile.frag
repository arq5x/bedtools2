BUILT_OBJECTS += obj/BamAncillary.o

obj/BamAncillary.o: src/utils/BamTools-Ancillary/BamAncillary.cpp obj/BamAncillary.d
	$(CXX_COMPILE)
