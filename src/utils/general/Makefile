OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/lineFileUtilities/ \
			-I$(UTILITIES_DIR)/version/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= ParseTools.h ParseTools.cpp PushBackStreamBuf.cpp PushBackStreamBuf.h CompressionTools.h CompressionTools.cpp \
		 Tokenizer.h Tokenizer.cpp CommonHelp.h CommonHelp.cpp ErrorMsg.h ErrorMsg.cpp
OBJECTS= ParseTools.o PushBackStreamBuf.o CompressionTools.o Tokenizer.o CommonHelp.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/string.o $(OBJ_DIR)/ParseTools.o $(OBJ_DIR)/PushBackStreamBuf.o $(OBJ_DIR)/Tokenizer.o $(OBJ_DIR)/CommonHelp.o

.PHONY: clean
