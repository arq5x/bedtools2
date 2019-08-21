BUILT_OBJECTS += obj/expand.o

obj/expand.o: src/expand/expand.cpp obj/expand.d
	$(CXX_COMPILE)
