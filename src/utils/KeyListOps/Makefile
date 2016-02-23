OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/general/ \
			-I$(UTILITIES_DIR)/lineFileUtilities/ \
			-I$(UTILITIES_DIR)/fileType/ \
			-I$(UTILITIES_DIR)/GenomeFile/ \
           -I$(UTILITIES_DIR)/FileRecordTools/ \
           -I$(UTILITIES_DIR)/FileRecordTools/FileReaders/ \
           -I$(UTILITIES_DIR)/FileRecordTools/Records/ \
           -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/src/ \
            -I$(UTILITIES_DIR)/version/
         
         
# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= KeyListOps.cpp KeyListOps.h KeyListOpsMethods.cpp KeyListOpsMethods.h
OBJECTS= KeyListOps.o KeyListOpsMethods.o
_EXT_OBJECTS=
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)


$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/KeyListOps.o $(OBJ_DIR)/KeyListOpsMethods.o

.PHONY: clean