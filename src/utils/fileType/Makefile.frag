BUILT_OBJECTS += obj/fileType.o obj/FileRecordTypeChecker.o

obj/fileType.o: src/utils/fileType/fileType.cpp obj/fileType.d
	$(CXX_COMPILE)

obj/FileRecordTypeChecker.o: src/utils/fileType/FileRecordTypeChecker.cpp obj/FileRecordTypeChecker.d
	$(CXX_COMPILE)
