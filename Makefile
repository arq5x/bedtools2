# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
# ==========================

# define our object and binary directories
export OBJ_DIR	= obj
export BIN_DIR	= bin
export SRC_DIR	= src
export CXX		= g++
export CXXFLAGS = -Wall -O2
export LIBS		= -lz



SUBDIRS = $(SRC_DIR)/annotateBed \
		  $(SRC_DIR)/bamToBed \
		  $(SRC_DIR)/bedToBam \
		  $(SRC_DIR)/bedToIgv \
		  $(SRC_DIR)/bed12ToBed6 \
		  $(SRC_DIR)/closestBed \
		  $(SRC_DIR)/complementBed \
		  $(SRC_DIR)/coverageBed \
		  $(SRC_DIR)/fastaFromBed \
		  $(SRC_DIR)/fjoin \
		  $(SRC_DIR)/genomeCoverageBed \
		  $(SRC_DIR)/intersectBed \
		  $(SRC_DIR)/linksBed \
		  $(SRC_DIR)/maskFastaFromBed \
		  $(SRC_DIR)/mergeBed	\
		  $(SRC_DIR)/overlap \
		  $(SRC_DIR)/pairToBed \
		  $(SRC_DIR)/pairToPair \
		  $(SRC_DIR)/shuffleBed \
		  $(SRC_DIR)/slopBed \
		  $(SRC_DIR)/sortBed \
		  $(SRC_DIR)/subtractBed \
		  $(SRC_DIR)/unionBedGraphs \
		  $(SRC_DIR)/windowBed

UTIL_SUBDIRS =	$(SRC_DIR)/utils/lineFileUtilities \
				$(SRC_DIR)/utils/bedFile \
				$(SRC_DIR)/utils/bedGraphFile \
				$(SRC_DIR)/utils/tabFile \
				$(SRC_DIR)/utils/genomeFile \
				$(SRC_DIR)/utils/gzstream \
				$(SRC_DIR)/utils/fileType \
				$(SRC_DIR)/utils/bedFilePE \
				$(SRC_DIR)/utils/sequenceUtilities \
				$(SRC_DIR)/utils/BamTools

all:
	[ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR)
	[ -d $(BIN_DIR) ] || mkdir -p $(BIN_DIR)
	
	@echo "Building BEDTools:"
	@echo "========================================================="
	
	@for dir in $(UTIL_SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done


.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
