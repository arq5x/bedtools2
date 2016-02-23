UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/
TOOL_DIR = ../../src/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/Contexts/ \
			-I$(UTILITIES_DIR)/general/ \
			-I$(UTILITIES_DIR)/fileType/ \
			-I$(UTILITIES_DIR)/lineFileUtilities/ \
           -I$(UTILITIES_DIR)/GenomeFile/ \
           -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/src \
           -I$(UTILITIES_DIR)/BlockedIntervals \
           -I$(UTILITIES_DIR)/BamTools-Ancillary \
           -I$(UTILITIES_DIR)/FileRecordTools/ \
           -I$(UTILITIES_DIR)/FileRecordTools/FileReaders/ \
           -I$(UTILITIES_DIR)/FileRecordTools/Records/ \
			-I$(UTILITIES_DIR)/KeyListOps/ \
           -I$(UTILITIES_DIR)/RecordOutputMgr/ \
           -I$(UTILITIES_DIR)/NewChromsweep \
           -I$(UTILITIES_DIR)/BinTree \
           -I$(UTILITIES_DIR)/version/ \
           -I$(UTILITIES_DIR)/ToolBase/ \
           -I$(TOOL_DIR)/intersectFile/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= subtractHelp.cpp subtractFile.cpp subtractFile.h
OBJECTS= subtractHelp.o subtractFile.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= subtractFile

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES)
	
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/subtractHelp.o $(OBJ_DIR)/subtractFile.o

.PHONY: clean
