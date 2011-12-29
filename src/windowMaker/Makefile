UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/genomeFile/ \
           -I$(UTILITIES_DIR)/BamTools/include

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= windowMakerMain.cpp windowMaker.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=genomeFile.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= windowMaker


all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/genomeFile/
	
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
