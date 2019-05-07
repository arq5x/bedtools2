BUILT_OBJECTS += obj/annotateMain.o obj/annotateBed.o

obj/annotateMain.o: src/annotateBed/annotateMain.cpp obj/annotateMain.d
	$(CXX_COMPILE)

obj/annotateBed.o: src/annotateBed/annotateBed.cpp obj/annotateBed.d
	$(CXX_COMPILE)
