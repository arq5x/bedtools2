BUILT_OBJECTS += obj/RecordOutputMgr.o

obj/RecordOutputMgr.o: src/utils/RecordOutputMgr/RecordOutputMgr.cpp obj/RecordOutputMgr.d
	$(CXX_COMPILE)
