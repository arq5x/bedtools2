UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/
TOOL_DIR = ../../src/

# -------------------
# define our includes
# -------------------

INCLUDES = -I$(UTILITIES_DIR)/bedFile/ \
			-I$(UTILITIES_DIR)/Contexts/ \
			-I$(UTILITIES_DIR)/general/ \
			-I$(UTILITIES_DIR)/fileType/ \
		    -I$(UTILITIES_DIR)/gzstream/ \
		   -I$(UTILITIES_DIR)/bedFile/ \
		   -I$(UTILITIES_DIR)/bedFilePE/ \
		   -I$(UTILITIES_DIR)/GenomeFile/ \
		   -I$(UTILITIES_DIR)/lineFileUtilities/ \
           -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/src \
           -I$(UTILITIES_DIR)/BlockedIntervals \
           -I$(UTILITIES_DIR)/BamTools-Ancillary \
           -I$(UTILITIES_DIR)/FileRecordTools/ \
           -I$(UTILITIES_DIR)/FileRecordTools/FileReaders/ \
           -I$(UTILITIES_DIR)/FileRecordTools/Records/ \
           -I$(UTILITIES_DIR)/RecordOutputMgr/ \
			-I$(UTILITIES_DIR)/KeyListOps/ \
           -I$(UTILITIES_DIR)/NewChromsweep \
           -I$(UTILITIES_DIR)/VectorOps \
 		   -I . \
           -I$(UTILITIES_DIR)/BinTree \
           -I$(UTILITIES_DIR)/version/ \
           -I$(UTILITIES_DIR)/ToolBase/ \
           -I$(TOOL_DIR)/intersectFile/ \
           -I$(TOOL_DIR)/jaccard/
            

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= fisherHelp.cpp fisher.cpp fisher.h kfunc.cpp \
		 ../utils/Contexts/ContextFisher.cpp \
		 ../utils/Contexts/ContextFisher.h
OBJECTS= fisherHelp.o fisher.o kfunc.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= Fisher

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(DFLAGS) $(INCLUDES)
	
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/fisherHelp.o $(OBJ_DIR)/fisher.o $(OBJ_DIR)/kfunc.o

.PHONY: clean
