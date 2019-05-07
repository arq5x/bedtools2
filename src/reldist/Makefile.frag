BUILT_OBJECTS += obj/reldistMain.o obj/reldist.o

obj/reldistMain.o: src/reldist/reldistMain.cpp obj/reldistMain.d
	$(CXX_COMPILE)

obj/reldist.o: src/reldist/reldist.cpp obj/reldist.d
	$(CXX_COMPILE)
