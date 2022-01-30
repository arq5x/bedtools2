
# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
# ==========================

SHELL := /bin/bash -e

VERSION_FILE=./src/utils/version/version_git.h
RELEASED_VERSION_FILE=./src/utils/version/version_release.txt


# define our object and binary directories
ifeq ($(VERBOSE),1)
CCPREFIX =
else
CCPREFIX = @
endif

OBJ_DIR	= obj
BIN_DIR	= bin
SRC_DIR	= src

CXX     = g++

PYTHON ?= $(shell python --version >/dev/null 2>&1 && echo "python" || echo python3)

ifeq ($(DEBUG),1)
BT_CPPFLAGS = -DDEBUG -D_DEBUG -D_FILE_OFFSET_BITS=64 -DWITH_HTS_CB_API $(INCLUDES)
BT_CXXFLAGS = -Wconversion -Wall -Wextra -g -O0
else
BT_CPPFLAGS = -D_FILE_OFFSET_BITS=64 -DWITH_HTS_CB_API $(INCLUDES)
BT_CXXFLAGS = -g -Wall -O2
endif

# If the user has specified to do so, tell the compile to use rand() (instead of mt19937).
ifeq ($(USE_RAND),1)
BT_CXXFLAGS += -DUSE_RAND
else
BT_CXXFLAGS += -std=c++11
endif

BT_LDFLAGS =
BT_LIBS    = -lz -lm -lbz2 -llzma -lpthread

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
		  $(SRC_DIR)/summaryFile \
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
				$(SRC_DIR)/utils/BamTools-Ancillary \
				$(SRC_DIR)/utils/BlockedIntervals \
				$(SRC_DIR)/utils/Fasta \
				$(SRC_DIR)/utils/VectorOps \
				$(SRC_DIR)/utils/GenomeFile \
				$(SRC_DIR)/utils/RecordOutputMgr \
				$(SRC_DIR)/utils/ToolBase \
				$(SRC_DIR)/utils/driver \
				$(SRC_DIR)/utils/BamTools

INCLUDES =		$(addprefix -I,$(SUBDIRS) $(UTIL_SUBDIRS)) \
				-I$(SRC_DIR)/utils/BamTools/include \
				-I$(HTSDIR) \
				-I$(SRC_DIR)/utils/lineFileUtilities \
				-I$(SRC_DIR)/utils/Point \
				-I$(SRC_DIR)/utils/stringUtilities


# Ensure that the user's $CXXFLAGS/etc are applied after bedtools's.
ALL_CXXFLAGS = $(BT_CXXFLAGS) $(CXXFLAGS)
ALL_CPPFLAGS = $(BT_CPPFLAGS) $(CPPFLAGS)
ALL_LDFLAGS  = $(BT_LDFLAGS) $(LDFLAGS)
ALL_LIBS     = $(BT_LIBS) $(LIBS)


all: print_banner $(BIN_DIR)/bedtools $(BIN_DIR)/intersectBed test/htsutil

static: print_banner $(BIN_DIR)/bedtools.static

BUILT_OBJECTS = $(OBJ_DIR)/bedtools.o
# Include all the Makefile fragments, which add to $(BUILT_OBJECTS)
include $(patsubst %,%/Makefile.frag,$(SUBDIRS) $(UTIL_SUBDIRS))

## Automatically generate C++ dependencies.
## $(DEPFLAGS) is a set of compiler flags that causes the compiler to generate
## dependencies as a byproduct (which we write to a temporary file, only moving
## it into place on successful compilations).
## We then include the dependency files into this Makefile.
## The subdirectories' Makefile fragments contain rules like the one for
## $(OBJ_DIR)/bedtools.o below, with a dependency on obj/foo.d so that if
## the dependency file is missing the target will be rebuilt.
##
DEPFLAGS = -MT $@ -MMD -MP -MF $*.Td

define CXX_COMPILE
@echo "  * compiling $<"
$(CCPREFIX) $(CC_WRAPPER) $(CXX) $(ALL_CXXFLAGS) $(ALL_CPPFLAGS) $(DEPFLAGS) -c -o $@ $<
@mv -f $*.Td $*.d
endef

$(OBJ_DIR)/%.d: ;
.PRECIOUS: $(OBJ_DIR)/%.d

-include $(patsubst %.o,%.d,$(BUILT_OBJECTS))

$(OBJ_DIR)/bedtools.o: $(SRC_DIR)/bedtools.cpp $(OBJ_DIR)/bedtools.d
	$(CXX_COMPILE)

$(OBJ_DIR)/htsutil.o: $(SRC_DIR)/htsutil.cpp $(OBJ_DIR)/htsutil.d
	$(CXX_COMPILE)

# HTSlib's htslib.mk provides rules to rebuild $(HTSDIR)/libhts.a.
HTSDIR = src/utils/htslib
include $(HTSDIR)/htslib.mk

# This order-only prerequisite ensures OBJ_DIR exists before building any .o file
# but ignores the directory's timestamp, which changes every time a .o file is written.
$(BUILT_OBJECTS): | $(OBJ_DIR)

$(BIN_DIR)/bedtools: autoversion $(BUILT_OBJECTS) $(HTSDIR)/libhts.a | $(BIN_DIR)
	@echo "- Building main bedtools binary."
	$(CCPREFIX) $(CC_WRAPPER) $(CXX) $(ALL_LDFLAGS) -o $(BIN_DIR)/bedtools $(BUILT_OBJECTS) $(HTSDIR)/libhts.a $(ALL_LIBS)
	@echo "done."

$(BIN_DIR)/bedtools.static: autoversion $(BUILT_OBJECTS) $(HTSDIR)/libhts.a | $(BIN_DIR)
	@echo "- Building main bedtools binary."
	$(CCPREFIX) $(CC_WRAPPER) $(CXX) -static $(ALL_LDFLAGS) -o $(BIN_DIR)/bedtools.static $(BUILT_OBJECTS) $(HTSDIR)/libhts.a $(ALL_LIBS)
	@echo "done."

test/htsutil: $(OBJ_DIR)/htsutil.o $(HTSDIR)/libhts.a
	$(CCPREFIX) $(CC_WRAPPER) $(CXX) $(ALL_LDFLAGS) -o test/htsutil $(OBJ_DIR)/htsutil.o $(HTSDIR)/libhts.a $(ALL_LIBS)

$(BIN_DIR)/intersectBed: | $(BIN_DIR)
	@echo "- Creating executables for old CLI."
	@$(PYTHON) scripts/makeBashScripts.py
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
	@echo "CXXFLAGS is [$(ALL_CXXFLAGS)]"
.PHONY: print_banner

# make the "obj/" and "bin/" directories, if they don't exist
$(OBJ_DIR) $(BIN_DIR):
	@mkdir -p $@


# Usually HTSlib's configure script has not been used (detected via config.mk
# not existing), so clean should also remove the generated config.h.
clean:
	@echo " * Cleaning up."
	@rm -f $(VERSION_FILE) $(OBJ_DIR)/* $(BIN_DIR)/* test/htsutil
	@cd src/utils/htslib && make clean > /dev/null
	@test -e src/utils/htslib/config.mk || rm -f src/utils/htslib/config.h
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
