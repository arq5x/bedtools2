BUILT_OBJECTS += obj/FileRecordMgr.o obj/FileRecordMergeMgr.o

obj/FileRecordMgr.o: src/utils/FileRecordTools/FileRecordMgr.cpp obj/FileRecordMgr.d
	$(CXX_COMPILE)

obj/FileRecordMergeMgr.o: src/utils/FileRecordTools/FileRecordMergeMgr.cpp obj/FileRecordMergeMgr.d
	$(CXX_COMPILE)
