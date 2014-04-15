OBJ_DIR = ../../../../obj/
BIN_DIR = ../../../../bin/
UTILITIES_DIR = ../../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/FileRecordTools/FileReaders/ \
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
SOURCES= RecordMgr.cpp RecordMgr.h Record.h Record.cpp EmptyRecord.h EmptyRecord.cpp Bed3Interval.h Bed3Interval.cpp \
	Bed4Interval.h Bed4Interval.cpp BedGraphInterval.h BedGraphInterval.cpp Bed5Interval.h Bed5Interval.cpp \
	Bed6Interval.h Bed6Interval.cpp \
	BedPlusInterval.h BedPlusInterval.cpp Bed12Interval.h Bed12Interval.cpp BamRecord.h BamRecord.cpp VcfRecord.h VcfRecord.cpp \
	GffRecord.h GffRecord.cpp RecordKeyList.h RecordKeyList.cpp BlockMgr.h BlockMgr.cpp StrandQueue.h StrandQueue.cpp
OBJECTS= RecordMgr.o Record.o EmptyRecord.o Bed3Interval.o Bed4Interval.o BedGraphInterval.o Bed5Interval.o Bed6Interval.o BedPlusInterval.o Bed12Interval.o BamRecord.o \
		VcfRecord.o GffRecord.o RecordKeyList.o BlockMgr.o StrandQueue.o
_EXT_OBJECTS=ParseTools.o QuickString.o ChromIdLookup.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/RecordMgr.o $(OBJ_DIR)/Record.o  $(OBJ_DIR)/EmptyRecord.o $(OBJ_DIR)/Bed3Interval.o $(OBJ_DIR)/Bed4Interval.o \
		$(OBJ_DIR)/BedGraphInterval.o $(OBJ_DIR)/Bed5Interval.o $(OBJ_DIR)/Bed6Interval.o \
		$(OBJ_DIR)/BedPlusInterval.o $(OBJ_DIR)/Bed12Interval.o $(OBJ_DIR)/BamRecord.o $(OBJ_DIR)/VcfRecord.o $(OBJ_DIR)/GffRecord.o $(OBJ_DIR)/BlockMgr.o $(OBJ_DIR)/StrandQueue.o

.PHONY: clean