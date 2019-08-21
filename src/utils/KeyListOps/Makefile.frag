BUILT_OBJECTS += obj/KeyListOps.o obj/KeyListOpsMethods.o

obj/KeyListOps.o: src/utils/KeyListOps/KeyListOps.cpp obj/KeyListOps.d
	$(CXX_COMPILE)

obj/KeyListOpsMethods.o: src/utils/KeyListOps/KeyListOpsMethods.cpp obj/KeyListOpsMethods.d
	$(CXX_COMPILE)
