BUILT_OBJECTS += obj/jaccardHelp.o obj/jaccard.o

obj/jaccardHelp.o: src/jaccard/jaccardHelp.cpp obj/jaccardHelp.d
	$(CXX_COMPILE)

obj/jaccard.o: src/jaccard/jaccard.cpp obj/jaccard.d
	$(CXX_COMPILE)
