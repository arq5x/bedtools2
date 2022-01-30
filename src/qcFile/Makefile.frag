BUILT_OBJECTS += obj/qcHelp.o obj/qcFile.o

obj/qcHelp.o: src/qcFile/qcHelp.cpp obj/qcHelp.d
	$(CXX_COMPILE)

obj/qcFile.o: src/qcFile/qcFile.cpp obj/qcFile.d
	$(CXX_COMPILE)
