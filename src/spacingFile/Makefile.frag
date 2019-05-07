BUILT_OBJECTS += obj/spacingHelp.o obj/spacingFile.o

obj/spacingHelp.o: src/spacingFile/spacingHelp.cpp obj/spacingHelp.d
	$(CXX_COMPILE)

obj/spacingFile.o: src/spacingFile/spacingFile.cpp obj/spacingFile.d
	$(CXX_COMPILE)
