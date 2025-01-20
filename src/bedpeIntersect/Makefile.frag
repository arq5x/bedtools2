BUILT_OBJECTS += obj/bedpeIntersect.o

obj/bedpeIntersect.o: src/bedpeIntersect/bedpeIntersect.cpp obj/bedpeIntersect.d
	$(CXX_COMPILE)
