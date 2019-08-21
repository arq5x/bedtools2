BUILT_OBJECTS += obj/Record.o obj/EmptyRecord.o obj/Bed3Interval.o obj/Bed4Interval.o obj/BedGraphInterval.o obj/Bed5Interval.o obj/Bed6Interval.o obj/PlusFields.o obj/GffRecord.o obj/GffPlusRecord.o obj/NoPosPlusRecord.o obj/BedPlusInterval.o obj/Bed12Interval.o obj/BamRecord.o obj/VcfRecord.o obj/BlockMgr.o obj/StrandQueue.o obj/RecordMgr.o obj/RecordList.o obj/RecordKeyList.o obj/RecordKeyVector.o

obj/Record.o: src/utils/FileRecordTools/Records/Record.cpp obj/Record.d
	$(CXX_COMPILE)

obj/EmptyRecord.o: src/utils/FileRecordTools/Records/EmptyRecord.cpp obj/EmptyRecord.d
	$(CXX_COMPILE)

obj/Bed3Interval.o: src/utils/FileRecordTools/Records/Bed3Interval.cpp obj/Bed3Interval.d
	$(CXX_COMPILE)

obj/Bed4Interval.o: src/utils/FileRecordTools/Records/Bed4Interval.cpp obj/Bed4Interval.d
	$(CXX_COMPILE)

obj/BedGraphInterval.o: src/utils/FileRecordTools/Records/BedGraphInterval.cpp obj/BedGraphInterval.d
	$(CXX_COMPILE)

obj/Bed5Interval.o: src/utils/FileRecordTools/Records/Bed5Interval.cpp obj/Bed5Interval.d
	$(CXX_COMPILE)

obj/Bed6Interval.o: src/utils/FileRecordTools/Records/Bed6Interval.cpp obj/Bed6Interval.d
	$(CXX_COMPILE)

obj/PlusFields.o: src/utils/FileRecordTools/Records/PlusFields.cpp obj/PlusFields.d
	$(CXX_COMPILE)

obj/GffRecord.o: src/utils/FileRecordTools/Records/GffRecord.cpp obj/GffRecord.d
	$(CXX_COMPILE)

obj/GffPlusRecord.o: src/utils/FileRecordTools/Records/GffPlusRecord.cpp obj/GffPlusRecord.d
	$(CXX_COMPILE)

obj/NoPosPlusRecord.o: src/utils/FileRecordTools/Records/NoPosPlusRecord.cpp obj/NoPosPlusRecord.d
	$(CXX_COMPILE)

obj/BedPlusInterval.o: src/utils/FileRecordTools/Records/BedPlusInterval.cpp obj/BedPlusInterval.d
	$(CXX_COMPILE)

obj/Bed12Interval.o: src/utils/FileRecordTools/Records/Bed12Interval.cpp obj/Bed12Interval.d
	$(CXX_COMPILE)

obj/BamRecord.o: src/utils/FileRecordTools/Records/BamRecord.cpp obj/BamRecord.d
	$(CXX_COMPILE)

obj/VcfRecord.o: src/utils/FileRecordTools/Records/VcfRecord.cpp obj/VcfRecord.d
	$(CXX_COMPILE)

obj/BlockMgr.o: src/utils/FileRecordTools/Records/BlockMgr.cpp obj/BlockMgr.d
	$(CXX_COMPILE)

obj/StrandQueue.o: src/utils/FileRecordTools/Records/StrandQueue.cpp obj/StrandQueue.d
	$(CXX_COMPILE)

obj/RecordMgr.o: src/utils/FileRecordTools/Records/RecordMgr.cpp obj/RecordMgr.d
	$(CXX_COMPILE)

obj/RecordList.o: src/utils/FileRecordTools/Records/RecordList.cpp obj/RecordList.d
	$(CXX_COMPILE)

obj/RecordKeyList.o: src/utils/FileRecordTools/Records/RecordKeyList.cpp obj/RecordKeyList.d
	$(CXX_COMPILE)

obj/RecordKeyVector.o: src/utils/FileRecordTools/Records/RecordKeyVector.cpp obj/RecordKeyVector.d
	$(CXX_COMPILE)
