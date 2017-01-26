OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include/ \
			-I$(UTILITIES_DIR)/general

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= fileType.cpp fileType.h FileRecordTypeChecker.cpp FileRecordTypeChecker.h
OBJECTS= fileType.o FileRecordTypeChecker.o
#_EXT_OBJECTS=ParseTools.o string.o
#EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/fileType.o $(OBJ_DIR)/FileRecordTypeChecker.o

.PHONY: clean