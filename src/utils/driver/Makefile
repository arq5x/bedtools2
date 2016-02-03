UTILITIES_DIR = ../../utils/
TOOLS_DIR = ../../../src/
OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/bedFile/ \
           -I$(UTILITIES_DIR)/lineFileUtilities/ \
			-I$(UTILITIES_DIR)/Contexts/ \
			-I$(UTILITIES_DIR)/general/ \
			-I$(UTILITIES_DIR)/fileType/ \
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
           -I$(UTILITIES_DIR)/version/ \
			-I$(TOOLS_DIR)/annotateBed/ \
			-I$(TOOLS_DIR)/bamToBed/ \
			-I$(TOOLS_DIR)/bamToFastq/ \
			-I$(TOOLS_DIR)/bedToBam/ \
			-I$(TOOLS_DIR)/bedpeToBam/ \
			-I$(TOOLS_DIR)/bedToIgv/ \
			-I$(TOOLS_DIR)/bed12ToBed6/ \
			-I$(TOOLS_DIR)/closestFile/ \
			-I$(TOOLS_DIR)/clusterBed/ \
			-I$(TOOLS_DIR)/complementFile/ \
			-I$(TOOLS_DIR)/coverageFile/ \
			-I$(TOOLS_DIR)/expand/ \
			-I$(TOOLS_DIR)/fastaFromBed/ \
			-I$(TOOLS_DIR)/flankBed/ \
			-I$(TOOLS_DIR)/genomeCoverageBed/ \
			-I$(TOOLS_DIR)/getOverlap/ \
			-I$(TOOLS_DIR)/groupBy/ \
			-I$(TOOLS_DIR)/intersectFile/ \
			-I$(TOOLS_DIR)/fisher/ \
			-I$(TOOLS_DIR)/jaccard/ \
			-I$(TOOLS_DIR)/linksBed/ \
			-I$(TOOLS_DIR)/maskFastaFromBed/ \
			-I$(TOOLS_DIR)/mapFile/ \
			-I$(TOOLS_DIR)/mergeFile/ \
			-I$(TOOLS_DIR)/multiBamCov/ \
			-I$(TOOLS_DIR)/multiIntersectBed/ \
			-I$(TOOLS_DIR)/nekSandbox1/ \
			-I$(TOOLS_DIR)/nucBed/ \
			-I$(TOOLS_DIR)/pairToBed/ \
			-I$(TOOLS_DIR)/pairToPair/ \
			-I$(TOOLS_DIR)/randomBed/ \
			-I$(TOOLS_DIR)/regressTest/ \
			-I$(TOOLS_DIR)/reldist/ \
			-I$(TOOLS_DIR)/sampleFile/ \
			-I$(TOOLS_DIR)/shiftBed/ \
			-I$(TOOLS_DIR)/shuffleBed/ \
			-I$(TOOLS_DIR)/slopBed/ \
			-I$(TOOLS_DIR)/sortBed/ \
			-I$(TOOLS_DIR)/spacingFile/ \
			-I$(TOOLS_DIR)/split/ \
			-I$(TOOLS_DIR)/subtractFile/ \
			-I$(TOOLS_DIR)/tagBam/ \
			-I$(TOOLS_DIR)/unionBedGraphs/ \
			-I$(TOOLS_DIR)/windowBed/ \
			-I$(TOOLS_DIR)/windowMaker/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= BedtoolsDriver.h BedtoolsDriver.cpp
OBJECTS= BedtoolsDriver.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES)
	
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/BedtoolsDriver.o

.PHONY: clean
