OBJ_DIR = ../../../../obj/
BIN_DIR = ../../../../bin/
UTILITIES_DIR = ../../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/lineFileUtilities/ \
		-I$(UTILITIES_DIR)/FileRecordTools/FileReaders/ \
		-I$(UTILITIES_DIR)/general/ \
		-I$(UTILITIES_DIR)/fileType/ \
		-I$(UTILITIES_DIR)/Contexts/ \
		-I$(UTILITIES_DIR)/GenomeFile/ \
		-I$(UTILITIES_DIR)/BamTools/include \
		-I$(UTILITIES_DIR)/BamTools/src/ \
		-I$(UTILITIES_DIR)/version/
			

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES = Record.cpp EmptyRecord.cpp Bed3Interval.cpp \
	Bed4Interval.cpp BedGraphInterval.cpp Bed5Interval.cpp \
	Bed6Interval.cpp PlusFields.cpp GffRecord.cpp GffPlusRecord.cpp \
	NoPosPlusRecord.cpp \
	BedPlusInterval.cpp Bed12Interval.cpp BamRecord.cpp VcfRecord.cpp \
	BlockMgr.cpp StrandQueue.cpp \
	RecordMgr.cpp RecordList.cpp RecordKeyList.cpp RecordKeyVector.cpp
OBJECTS= Record.o EmptyRecord.o Bed3Interval.o Bed4Interval.o BedGraphInterval.o Bed5Interval.o Bed6Interval.o PlusFields.o BedPlusInterval.o Bed12Interval.o BamRecord.o \
	GffRecord.o GffPlusRecord.o VcfRecord.o  NoPosPlusRecord.o BlockMgr.o StrandQueue.o RecordMgr.o RecordList.o RecordKeyList.o RecordKeyVector.o
_EXT_OBJECTS=ParseTools.o string.o ChromIdLookup.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/RecordMgr.o $(OBJ_DIR)/RecordList.o $(OBJ_DIR)/Record.o  $(OBJ_DIR)/EmptyRecord.o $(OBJ_DIR)/Bed3Interval.o $(OBJ_DIR)/Bed4Interval.o \
		$(OBJ_DIR)/BedGraphInterval.o $(OBJ_DIR)/Bed5Interval.o $(OBJ_DIR)/Bed6Interval.o $(OBJ_DIR)/PlusFields.o $(OBJ_DIR)/GffPlusRecord.o $(OBJ_DIR)/NoPosPlusRecord.o \
		$(OBJ_DIR)/BedPlusInterval.o $(OBJ_DIR)/Bed12Interval.o $(OBJ_DIR)/BamRecord.o $(OBJ_DIR)/VcfRecord.o $(OBJ_DIR)/GffRecord.o $(OBJ_DIR)/BlockMgr.o $(OBJ_DIR)/StrandQueue.o \
		$(OBJ_DIR)/RecordKeyList.o $(OBJ_DIR)/RecordKeyVector.o

.PHONY: clean