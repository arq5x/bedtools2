BUILT_OBJECTS += obj/BinTree.o

obj/BinTree.o: src/utils/BinTree/BinTree.cpp obj/BinTree.d
	$(CXX_COMPILE)
