OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/general/ \
		 	-I$(UTILITIES_DIR)/lineFileUtilities/ \
 			-I$(UTILITIES_DIR)/KeyListOps/ \
           -I$(UTILITIES_DIR)/BamTools/include/	

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= NewGenomeFile.cpp NewGenomeFile.h GenomeFile.cpp GenomeFile.h
OBJECTS= NewGenomeFile.o GenomeFile.o
#_EXT_OBJECTS=lineFileUtilities.o fileType.o
#EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

$(BUILT_OBJECTS): $(SOURCES) $(SUBDIRS)
	@echo "  * compiling GenomeFile.cpp"
	@$(CXX) -c -o $(OBJ_DIR)/GenomeFile.o GenomeFile.cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)
	@echo "  * compiling NewGenomeFile.cpp"
	@$(CXX) -c -o $(OBJ_DIR)/NewGenomeFile.o NewGenomeFile.cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)

#$(EXT_OBJECTS):
#	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/lineFileUtilities/
#	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/fileType/

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/NewGenomeFile.o $(OBJ_DIR)/GenomeFile.o

.PHONY: clean
