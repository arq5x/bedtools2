BUILT_OBJECTS += obj/NewGenomeFile.o obj/GenomeFile.o

obj/NewGenomeFile.o: src/utils/GenomeFile/NewGenomeFile.cpp obj/NewGenomeFile.d
	$(CXX_COMPILE)

obj/GenomeFile.o: src/utils/GenomeFile/GenomeFile.cpp obj/GenomeFile.d
	$(CXX_COMPILE)
