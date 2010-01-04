# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
# ==========================

# define our object and binary directories
export OBJ_DIR = obj
export BIN_DIR = bin
export SRC_DIR = src

# define some default flags
export CFLAGS ?= -Wall -O3
export CXXFLAGS ?= $(CFLAGS)
export LDFLAGS ?= -Wl,-s
export CXX ?= g++

# define our platform
#export BLD_PLATFORM ?= linux64-core2
#include includes/$(BLD_PLATFORM).inc

# define our source subdirectories
SUBDIRS = $(SRC_DIR)/closestBed $(SRC_DIR)/complementBed $(SRC_DIR)/coverageBed $(SRC_DIR)/intersectBed $(SRC_DIR)/mergeBed $(SRC_DIR)/genomeCoverageBed $(SRC_DIR)/fastaFromBed $(SRC_DIR)/shuffleBed $(SRC_DIR)/slopBed $(SRC_DIR)/sortBed $(SRC_DIR)/windowBed $(SRC_DIR)/subtractBed $(SRC_DIR)/linksBed $(SRC_DIR)/pairToBed $(SRC_DIR)/pairToPair $(SRC_DIR)/maskFastaFromBed $(SRC_DIR)/bamToBed

UTIL_SUBDIRS =  $(SRC_DIR)/utils/lineFileUtilities $(SRC_DIR)/utils/bedFile $(SRC_DIR)/utils/bedFilePE $(SRC_DIR)/utils/sequenceUtilities

all:

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
