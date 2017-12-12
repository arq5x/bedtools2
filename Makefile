
# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
# ==========================

SHELL := /bin/bash -e

VERSION_FILE=./src/utils/version/version_git.h
RELEASED_VERSION_FILE=./src/utils/version/version_release.txt



# define our object and binary directories
export OBJ_DIR	= obj
export BIN_DIR	= bin
export SRC_DIR	= src
export UTIL_DIR	= src/utils
export CXX		= g++
ifeq ($(DEBUG),1)
export CXXFLAGS = -Wall -Wextra -DDEBUG -D_DEBUG -g -O0 -D_FILE_OFFSET_BITS=64 -fPIC $(INCLUDES)
else
export CXXFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC $(INCLUDES)
endif

# If the user has specified to do so, tell the compile to use rand() (instead of mt19937).
ifeq ($(USE_RAND),1)
export CXXFLAGS += -DUSE_RAND
else
export CXXFLAGS += -std=c++11
endif 

export LIBS		= -lz
export BT_ROOT  = src/utils/BamTools/

prefix ?= /usr/local

SUBDIRS = $(SRC_DIR)/annotateBed \
		  $(SRC_DIR)/bamToBed \
		  $(SRC_DIR)/bamToFastq \
		  $(SRC_DIR)/bedToBam \
		  $(SRC_DIR)/bedpeToBam \
		  $(SRC_DIR)/bedToIgv \
		  $(SRC_DIR)/bed12ToBed6 \
		  $(SRC_DIR)/closestFile \
		  $(SRC_DIR)/clusterBed \
		  $(SRC_DIR)/complementFile \
		  $(SRC_DIR)/coverageFile \
		  $(SRC_DIR)/expand \
		  $(SRC_DIR)/fastaFromBed \
		  $(SRC_DIR)/flankBed \
		  $(SRC_DIR)/genomeCoverageBed \
		  $(SRC_DIR)/getOverlap \
		  $(SRC_DIR)/groupBy \
		  $(SRC_DIR)/intersectFile \
		  $(SRC_DIR)/fisher \
		  $(SRC_DIR)/jaccard \
		  $(SRC_DIR)/linksBed \
		  $(SRC_DIR)/maskFastaFromBed \
		  $(SRC_DIR)/mapFile \
		  $(SRC_DIR)/mergeFile \
		  $(SRC_DIR)/multiBamCov \
		  $(SRC_DIR)/multiIntersectBed \
		  $(SRC_DIR)/nucBed \
		  $(SRC_DIR)/pairToBed \
		  $(SRC_DIR)/pairToPair \
		  $(SRC_DIR)/randomBed \
		  $(SRC_DIR)/regressTest \
		  $(SRC_DIR)/reldist \
		  $(SRC_DIR)/sampleFile \
		  $(SRC_DIR)/shiftBed \
		  $(SRC_DIR)/shuffleBed \
		  $(SRC_DIR)/slopBed \
		  $(SRC_DIR)/sortBed \
		  $(SRC_DIR)/spacingFile \
		  $(SRC_DIR)/split \
		  $(SRC_DIR)/subtractFile \
		  $(SRC_DIR)/tagBam \
		  $(SRC_DIR)/unionBedGraphs \
		  $(SRC_DIR)/windowBed \
		  $(SRC_DIR)/windowMaker

UTIL_SUBDIRS =	$(SRC_DIR)/utils/FileRecordTools \
				$(SRC_DIR)/utils/FileRecordTools/FileReaders \
				$(SRC_DIR)/utils/FileRecordTools/Records \
				$(SRC_DIR)/utils/bedFile \
				$(SRC_DIR)/utils/BinTree \
				$(SRC_DIR)/utils/version \
				$(SRC_DIR)/utils/bedGraphFile \
				$(SRC_DIR)/utils/chromsweep \
				$(SRC_DIR)/utils/Contexts \
				$(SRC_DIR)/utils/general \
				$(SRC_DIR)/utils/gzstream \
				$(SRC_DIR)/utils/fileType \
				$(SRC_DIR)/utils/bedFilePE \
				$(SRC_DIR)/utils/KeyListOps \
				$(SRC_DIR)/utils/NewChromsweep \
				$(SRC_DIR)/utils/sequenceUtilities \
				$(SRC_DIR)/utils/tabFile \
				$(SRC_DIR)/utils/BamTools \
				$(SRC_DIR)/utils/BamTools-Ancillary \
				$(SRC_DIR)/utils/BlockedIntervals \
				$(SRC_DIR)/utils/Fasta \
				$(SRC_DIR)/utils/VectorOps \
				$(SRC_DIR)/utils/GenomeFile \
				$(SRC_DIR)/utils/RecordOutputMgr \
				$(SRC_DIR)/utils/ToolBase \
				$(SRC_DIR)/utils/driver

BUILT_OBJECTS = $(OBJ_DIR)/*.o


INCLUDES =	-I$(SRC_DIR)/utils/bedFile \
				-I$(SRC_DIR)/utils/BinTree \
				-I$(SRC_DIR)/utils/version \
				-I$(SRC_DIR)/utils/bedGraphFile \
				-I$(SRC_DIR)/utils/chromsweep \
				-I$(SRC_DIR)/utils/Contexts \
				-I$(SRC_DIR)/utils/FileRecordTools \
				-I$(SRC_DIR)/utils/FileRecordTools/FileReaders \
				-I$(SRC_DIR)/utils/FileRecordTools/Records \
				-I$(SRC_DIR)/utils/general \
				-I$(SRC_DIR)/utils/gzstream \
				-I$(SRC_DIR)/utils/fileType \
				-I$(SRC_DIR)/utils/gzstream/ \
				-I$(SRC_DIR)/utils/lineFileUtilities \
				-I$(SRC_DIR)/utils/KeyListOps \
				-I$(SRC_DIR)/utils/NewChromsweep \
				-I$(SRC_DIR)/utils/sequenceUtilities \
				-I$(SRC_DIR)/utils/tabFile \
				-I$(SRC_DIR)/utils/BamTools \
				-I$(SRC_DIR)/utils/BamTools/include \
				-I$(SRC_DIR)/utils/BamTools/src \
				-I$(SRC_DIR)/utils/BamTools-Ancillary \
				-I$(SRC_DIR)/utils/BlockedIntervals \
				-I$(SRC_DIR)/utils/Fasta \
				-I$(SRC_DIR)/utils/VectorOps \
				-I$(SRC_DIR)/utils/GenomeFile \
				-I$(SRC_DIR)/utils/RecordOutputMgr \
				-I$(SRC_DIR)/utils/ToolBase \
				-I$(SRC_DIR)/utils/driver \


all: print_banner $(OBJ_DIR) $(BIN_DIR) autoversion $(UTIL_SUBDIRS) $(SUBDIRS)
	@echo "- Building main bedtools binary."
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c src/bedtools.cpp -o obj/bedtools.o $(INCLUDES)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(BIN_DIR)/bedtools $(BUILT_OBJECTS) -L$(UTIL_DIR)/BamTools/lib/ -lbamtools $(LIBS) $(LDFLAGS) $(INCLUDES)
	@echo "done."

	@echo "- Creating executables for old CLI."
	@python scripts/makeBashScripts.py
	@chmod +x bin/*
	@echo "done."


.PHONY: all

install: all
	mkdir -p $(DESTDIR)$(prefix)/bin
	for file in bin/* ; do \
		cp -f $$file $(DESTDIR)$(prefix)/bin; \
	done

print_banner:
	@echo "Building BEDTools:"
	@echo "========================================================="
	$(info $$CXXFLAGS is [${CXXFLAGS}])
.PHONY: print_banner

# make the "obj/" and "bin/" directories, if they don't exist
$(OBJ_DIR) $(BIN_DIR):
	@mkdir -p $@


# One special case: All (or almost all) programs requires the BamTools API files to be created first.
.PHONY: bamtools_api
bamtools_api:
	@$(MAKE) --no-print-directory --directory=$(BT_ROOT) api
$(UTIL_SUBDIRS) $(SUBDIRS): bamtools_api


# even though these are real directories, treat them as phony targets, forcing to always go in them are re-make.
# a future improvement would be the check for the compiled object, and rebuild only if the source code is newer.
.PHONY: $(UTIL_SUBDIRS) $(SUBDIRS)
$(UTIL_SUBDIRS) $(SUBDIRS): $(OBJ_DIR) $(BIN_DIR)
	@echo "- Building in $@"
	@$(MAKE) --no-print-directory --directory=$@

clean:
	@$(MAKE) --no-print-directory --directory=$(BT_ROOT) clean_api
	@echo " * Cleaning up."
	@rm -f $(VERSION_FILE) $(OBJ_DIR)/* $(BIN_DIR)/*
.PHONY: clean

test: all
	@cd test; bash test.sh

.PHONY: test


## For BEDTools developers (not users):
## When you want to release (and tag) a new version, run:
##   $ make setversion VERSION=v2.17.2
## This will:
##   1. Update the "/src/utils/version/version_release.txt" file
##   2. Commit the file
##   3. Git-Tag the commit with the latest version
##
.PHONY: setversion
setversion:
    ifeq "$(VERSION)" ""
		$(error please set VERSION variable to the new version (e.g "make setversion VERSION=v2.17.2"))
    endif
		@echo "# This file was auto-generated by running \"make setversion VERSION=$(VERSION)\"" > "$(RELEASED_VERSION_FILE)"
		@echo "# on $$(date) ." >> "$(RELEASED_VERSION_FILE)"
		@echo "# Please do not edit or commit this file manually." >> "$(RELEASED_VERSION_FILE)"
		@echo "#" >> "$(RELEASED_VERSION_FILE)"
		@echo "$(VERSION)" >> $(RELEASED_VERSION_FILE)
		@git add $(RELEASED_VERSION_FILE)
		@git commit -q -m "Setting Release-Version $(VERSION)"
		@git tag "$(VERSION)"
		@echo "Version updated to $(VERSION)."
		@echo ""
		@echo "Don't forget to push the commits AND the tags:"
		@echo "  git push --all --tags"
		@echo ""


## Automatic version detection
##
## What's going on here?
## 1. If there's a ".git" repository - use the version from the repository.
##    ignore any released-version file. git repository is authorative.
##
## 2, If there's no ".git" repository,
##    get the "released" version number from the release-version file.
##
##    2.1. If the current directory looks like "arq5x-bedtools-XXXXXXX",
##         assume "-XXXXXX" is the last revision number (and the user
##         probably downloaded the ZIP from github).
##         Append the revision number to the released version string.
##
## 3. Compare the detected version (from steps 1,2) to the current string
##    in ./src/utils/version/version_git.h .
##    If they differ, update the header file - will cause a recompilation
##    of version.o .
##
.PHONY: autoversion
autoversion:
	@( \
	if [ -d ".git" ] && which git > /dev/null ; then \
		DETECTED_VERSION=$$(git describe --always --tags --dirty) ; \
	else \
		DETECTED_VERSION=$$(grep -v "^#" "$(RELEASED_VERSION_FILE)") ; \
		if basename $$(pwd) | grep -q "^[[:alnum:]]*-bedtools-[[:alnum:]]*$$" ; then \
			DETECTED_VERSION=$${DETECTED_VERSION}-zip-$$(basename "$$(pwd)" | sed 's/^[[:alnum:]]*-bedtools-//') ; \
		fi ; \
	fi ; \
	\
	CURRENT_VERSION="" ; \
	[ -e "$(VERSION_FILE)" ] && CURRENT_VERSION=$$(grep "define VERSION_GIT " "$(VERSION_FILE)" | cut -f3 -d" " | sed 's/"//g') ; \
	\
	echo "DETECTED_VERSION = $$DETECTED_VERSION" ; \
	echo "CURRENT_VERSION  = $$CURRENT_VERSION" ; \
	if [ "$${DETECTED_VERSION}" != "$${CURRENT_VERSION}" ] ; then \
		echo "Updating version file." ; \
		echo "#ifndef VERSION_GIT_H" > $(VERSION_FILE) ; \
		echo "#define VERSION_GIT_H" >> $(VERSION_FILE) ; \
		echo "#define VERSION_GIT \"$${DETECTED_VERSION}\"" >> $(VERSION_FILE) ; \
		echo "#endif /* VERSION_GIT_H */" >> $(VERSION_FILE) ; \
	fi )
