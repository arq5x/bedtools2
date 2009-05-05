CXX = g++
CXXFLAGS = -O3  -Wall
LDFLAGS = 
LIBFLAGS = -lc

UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/bedFile/ -I$(UTILITIES_DIR)/lineFileUtilities/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= complementMain.cpp complementBed.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=bedFile.o lineFileUtilities.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= complementBed
LIBS=

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LIBFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bedFile/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/lineFileUtilities/
	
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean