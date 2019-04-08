BUILT_OBJECTS += obj/fisherHelp.o obj/fisher.o obj/kfunc.o

obj/fisherHelp.o: src/fisher/fisherHelp.cpp obj/fisherHelp.d
	$(CXX_COMPILE)

obj/fisher.o: src/fisher/fisher.cpp obj/fisher.d
	$(CXX_COMPILE)

obj/kfunc.o: src/fisher/kfunc.cpp obj/kfunc.d
	$(CXX_COMPILE)
