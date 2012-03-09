UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include \
		   -I$(UTILITIES_DIR)/sequenceUtilities/ \
		   -I$(UTILITIES_DIR)/version/


# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= bamToFastqMain.cpp bamToFastq.cpp bamToFastq.h
OBJECTS= bamToFastqMain.o bamToFastq.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))


all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) 

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/bamToFastqMain.o  $(OBJ_DIR)/bamToFastq.o

.PHONY: clean
