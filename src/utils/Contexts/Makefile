OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES =  -I$(UTILITIES_DIR)/general/ \
			-I$(UTILITIES_DIR)/fileType/ \
			-I$(UTILITIES_DIR)/FileRecordTools/FileReaders/ \
			-I$(UTILITIES_DIR)/GenomeFile/ \
			-I$(UTILITIES_DIR)/BamTools/include/ \
			-I$(UTILITIES_DIR)/version/
			

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= Context.cpp Context.h
OBJECTS= Context.o
_EXT_OBJECTS=ParseTools.o QuickString.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/Context.o

.PHONY: clean