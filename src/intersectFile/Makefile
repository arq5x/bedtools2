UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/Contexts/ \
			-I$(UTILITIES_DIR)/general/ \
			-I$(UTILITIES_DIR)/fileType/ \
			-I$(UTILITIES_DIR)/lineFileUtilities/ \
			-I$(UTILITIES_DIR)/gzstream/ \
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
           -I$(UTILITIES_DIR)/ToolBase \
           -I$(UTILITIES_DIR)/version/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= intersectHelp.cpp intersectFile.cpp intersectFile.h
OBJECTS= intersectHelp.o intersectFile.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= intersectFile

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(DFLAGS) $(INCLUDES)
	
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/intersectHelp.o $(OBJ_DIR)/intersectFile.o

.PHONY: clean
