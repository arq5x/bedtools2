BUILT_OBJECTS += obj/FileReader.o obj/SingleLineDelimTextFileReader.o obj/BamFileReader.o obj/BufferedStreamMgr.o obj/InputStreamMgr.o

obj/FileReader.o: src/utils/FileRecordTools/FileReaders/FileReader.cpp obj/FileReader.d
	$(CXX_COMPILE)

obj/SingleLineDelimTextFileReader.o: src/utils/FileRecordTools/FileReaders/SingleLineDelimTextFileReader.cpp obj/SingleLineDelimTextFileReader.d
	$(CXX_COMPILE)

obj/BamFileReader.o: src/utils/FileRecordTools/FileReaders/BamFileReader.cpp obj/BamFileReader.d
	$(CXX_COMPILE)

obj/BufferedStreamMgr.o: src/utils/FileRecordTools/FileReaders/BufferedStreamMgr.cpp obj/BufferedStreamMgr.d
	$(CXX_COMPILE)

obj/InputStreamMgr.o: src/utils/FileRecordTools/FileReaders/InputStreamMgr.cpp obj/InputStreamMgr.d
	$(CXX_COMPILE)
