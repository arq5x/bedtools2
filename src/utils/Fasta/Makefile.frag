BUILT_OBJECTS += obj/Fasta.o obj/split.o

obj/Fasta.o: src/utils/Fasta/Fasta.cpp obj/Fasta.d
	$(CXX_COMPILE)

obj/split.o: src/utils/Fasta/split.cpp obj/split.d
	$(CXX_COMPILE)
