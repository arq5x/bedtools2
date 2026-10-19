
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
BT_LDFLAGS =
else
BT_CPPFLAGS = -D_FILE_OFFSET_BITS=64 -DWITH_HTS_CB_API $(INCLUDES)
ifeq ($(FAST),1)
# The optimized build, `make FAST=1`: the same code built with -O3, link-time
# optimization (LTO) and profile-guided optimization (PGO), each of LTO and PGO
# used when the compiler passes its probe below and skipped when it does not.
# Its output is the same as the default build's. A plain `make` builds with
# -O2 and neither. Run `make clean` when switching between the two, and give
# FAST=1 to every make on a tree built with it (`make FAST=1 test`, `static`,
# `install`), since each of them links bin/bedtools again. DEBUG=1 ignores it.
BT_CXXFLAGS = -g -Wall -O3
BT_LDFLAGS = -g -O3
# Link-time optimization inlines calls across source files. It is used when the
# compiler and linker accept it (clang on Linux also needs lld or the LLVMgold
# plugin); `make FAST=1 LTO=` builds without it.
ifeq ($(origin LTO),undefined)
LTO := $(shell d=`mktemp -d 2>/dev/null || mktemp -d -t bedtools` && \
	printf 'int main(void) { return 0; }\n' > $$d/t.cpp && \
	for f in -flto=auto -flto; do \
		if $(CXX) $(CXXFLAGS) $$f $(LDFLAGS) -o $$d/t $$d/t.cpp > /dev/null 2>&1; then echo $$f; break; fi; \
	done; rm -f $$d/t $$d/t.cpp; rmdir $$d)
endif
BT_CXXFLAGS += $(LTO)
BT_LDFLAGS += $(LTO)
# Profile-guided optimization: `make FAST=1` compiles an instrumented copy of
# every object into obj/pgo/, runs scripts/pgo-train.sh with it on inputs from
# test/, and compiles the release objects with the profile it records.
# Functions the training never runs are optimized as they would be without a
# profile. It is used when the compiler can build, run and read back a program
# instrumented this way (GCC; clang writes .profraw files, which need
# llvm-profdata, and falls back); `make FAST=1 PGO=` builds without it.
ifeq ($(origin PGO),undefined)
PGO := $(shell d=`mktemp -d 2>/dev/null || mktemp -d -t bedtools` && { \
	printf 'int f();\nstatic int g = f();\nint f() { return 1; }\nint main(void) { return g - 1; }\n' > $$d/t.cpp && \
	$(CXX) $(CXXFLAGS) $(LTO) -fprofile-generate -dumpbase $$d/u -c -o $$d/t.o $$d/t.cpp > /dev/null 2>&1 && \
	$(CXX) $(CXXFLAGS) $(LTO) -fprofile-generate $(LDFLAGS) -o $$d/t $$d/t.o > /dev/null 2>&1 && \
	(cd $$d && ./t) > /dev/null 2>&1 && test -f $$d/u.gcda && \
	$(CXX) $(CXXFLAGS) $(LTO) -fprofile-use -fprofile-partial-training -Werror -c -o $$d/u.o $$d/t.cpp > /dev/null 2>&1 && \
	$(CXX) $(CXXFLAGS) $(LTO) -fprofile-use -fprofile-partial-training $(LDFLAGS) -o $$d/u $$d/u.o > /dev/null 2>&1 && \
	echo 1; rm -f $$d/t $$d/t.o $$d/t.cpp $$d/u $$d/u.o $$d/u.gcda $$d/*.profraw; rmdir $$d; })
endif
ifneq ($(PGO),)
PGO_GEN_FLAGS = -fprofile-generate
PGO_USE_FLAGS = -fprofile-use -fprofile-partial-training
endif
else
BT_CXXFLAGS = -g -Wall -O2
BT_LDFLAGS =
endif
endif

# If the user has specified to do so, tell the compile to use rand() (instead of mt19937).
ifeq ($(USE_RAND),1)
BT_CXXFLAGS += -DUSE_RAND
else
BT_CXXFLAGS += -std=c++11
endif

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

ifneq ($(PGO_USE_FLAGS),)
# Profile-guided optimization (see PGO above). The instrumented objects are
# compiled from the same sources as $(BUILT_OBJECTS), found through vpath, and
# with the same auxiliary name (-dumpbase), which GCC uses to tell apart the
# profiles of functions local to an object and to name its .gcda file. The
# release objects and binaries take the profile, after bedtools's flags and
# before the user's.
ALL_CXXFLAGS = $(BT_CXXFLAGS) $(PGO_USE_FLAGS) $(CXXFLAGS)
ALL_LDFLAGS  = $(BT_LDFLAGS) $(PGO_USE_FLAGS) $(LDFLAGS)
PGO_DIR     = $(OBJ_DIR)/pgo
PGO_OBJECTS = $(patsubst $(OBJ_DIR)/%,$(PGO_DIR)/%,$(BUILT_OBJECTS))
vpath %.cpp $(SRC_DIR) $(SUBDIRS) $(UTIL_SUBDIRS)
vpath %.C $(SRC_DIR)/utils/gzstream

define PGO_COMPILE
@echo "  * compiling $< (instrumented)"
$(CCPREFIX) $(CC_WRAPPER) $(CXX) $(BT_CXXFLAGS) $(PGO_GEN_FLAGS) $(CXXFLAGS) $(ALL_CPPFLAGS) -dumpbase $(OBJ_DIR)/$* -MT $@ -MMD -MP -MF $(PGO_DIR)/$*.Td -c -o $@ $<
@mv -f $(PGO_DIR)/$*.Td $(PGO_DIR)/$*.d
endef

$(PGO_DIR)/%.o: %.cpp $(PGO_DIR)/%.d | $(PGO_DIR) autoversion
	$(PGO_COMPILE)

$(PGO_DIR)/%.o: %.C $(PGO_DIR)/%.d | $(PGO_DIR) autoversion
	$(PGO_COMPILE)

-include $(patsubst %.o,%.d,$(PGO_OBJECTS))

$(PGO_DIR):
	@mkdir -p $@

$(PGO_DIR)/bedtools: $(PGO_OBJECTS) $(HTSDIR)/libhts.a
	@echo "- Building the instrumented bedtools binary."
	$(CCPREFIX) $(CC_WRAPPER) $(CXX) $(BT_LDFLAGS) $(PGO_GEN_FLAGS) $(LDFLAGS) -o $@ $(PGO_OBJECTS) $(HTSDIR)/libhts.a $(ALL_LIBS)

# Each training call adds its counts to one .gcda file per object, written
# under GCOV_PREFIX. A release object's profile, $(OBJ_DIR)/<object>.gcda, is
# replaced only when the new one differs, so only those objects recompile.
$(PGO_DIR)/profile: $(PGO_DIR)/bedtools scripts/pgo-train.sh
	@echo "- Training the instrumented binary (scripts/pgo-train.sh)."
	@rm -rf $(PGO_DIR)/train
	$(CCPREFIX) GCOV_PREFIX=$(CURDIR)/$(PGO_DIR)/train/gcda sh scripts/pgo-train.sh $(PGO_DIR)/bedtools $(PGO_DIR)/train/data || \
		{ echo "The training failed; \`make FAST=1 PGO=\` builds without profile-guided optimization." >&2; exit 1; }
	@for f in `find $(PGO_DIR)/train/gcda -name '*.gcda'`; do \
		o=$(OBJ_DIR)/`basename $$f`; cmp -s $$f $$o || cp $$f $$o; \
	done
	@rm -rf $(PGO_DIR)/train
	@touch $@

# Release objects compile after the training; one whose profile changed
# recompiles. htsutil (test/htsutil) is not trained.
$(BUILT_OBJECTS): | $(PGO_DIR)/profile
PGO_GCDA := $(wildcard $(patsubst %.o,%.gcda,$(BUILT_OBJECTS)))
ifneq ($(PGO_GCDA),)
$(PGO_GCDA): $(PGO_DIR)/profile ;
$(foreach g,$(PGO_GCDA),$(eval $(patsubst %.gcda,%.o,$(g)): $(g)))
endif
$(OBJ_DIR)/htsutil.o test/htsutil: PGO_USE_FLAGS =
endif

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
	@rm -rf $(OBJ_DIR)/pgo
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
