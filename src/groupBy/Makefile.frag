BUILT_OBJECTS += obj/groupByHelp.o obj/groupBy.o

obj/groupByHelp.o: src/groupBy/groupByHelp.cpp obj/groupByHelp.d
	$(CXX_COMPILE)

obj/groupBy.o: src/groupBy/groupBy.cpp obj/groupBy.d
	$(CXX_COMPILE)
