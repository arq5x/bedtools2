OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../

INCLUDES = -I$(UTILITIES_DIR)/bedFile/ -I$(UTILITIES_DIR)/lineFileUtilities/ -I$(UTILITIES_DIR)/gzstream/ -I$(UTILITIES_DIR)/fileType/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= BamReader.cpp BamWriter.cpp	BGZF.cpp BamAncillary.cpp BamIndex.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean