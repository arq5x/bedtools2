OBJ_DIR = ../../../../obj/
BIN_DIR = ../../../../bin/
UTILITIES_DIR = ../../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include/ \
			-I$(UTILITIES_DIR)/BamTools/src/ \
			-I$(UTILITIES_DIR)/fileType/ \
			-I$(UTILITIES_DIR)/Contexts/ \
			-I$(UTILITIES_DIR)/GenomeFile/ \
			-I$(UTILITIES_DIR)/general/ \
			-I$(UTILITIES_DIR)/gzstream/ \
            -I$(UTILITIES_DIR)/version/ \
            -I$(UTILITIES_DIR)/lineFileUtilities/ \
			

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= FileReader.h FileReader.cpp  \
		SingleLineDelimTextFileReader.h SingleLineDelimTextFileReader.cpp BamFileReader.h BamFileReader.cpp \
		BufferedStreamMgr.h BufferedStreamMgr.cpp InputStreamMgr.h InputStreamMgr.cpp
OBJECTS= FileReader.o SingleLineDelimTextFileReader.o BamFileReader.o BufferedStreamMgr.o InputStreamMgr.o 
_EXT_OBJECTS=string.o CompressionTools.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/FileReader.o $(OBJ_DIR)/SingleLineDelimTextFileReader.o \
			$(OBJ_DIR)/BinaryFileReader.o $(OBJ_DIR)/BamFileReader.o $(OBJ_DIR)/BufferedStreamMgr.o $(OBJ_DIR)/InputStreamMgr.o $(OBJ_DIR)/PushBackGzStream.o

.PHONY: clean
