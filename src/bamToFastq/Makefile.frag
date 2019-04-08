BUILT_OBJECTS += obj/bamToFastqMain.o obj/bamToFastq.o

obj/bamToFastqMain.o: src/bamToFastq/bamToFastqMain.cpp obj/bamToFastqMain.d
	$(CXX_COMPILE)

obj/bamToFastq.o: src/bamToFastq/bamToFastq.cpp obj/bamToFastq.d
	$(CXX_COMPILE)
