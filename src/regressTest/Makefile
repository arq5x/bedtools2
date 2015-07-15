UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/bedFile/ \
			-I$(UTILITIES_DIR)/Contexts/ \
           -I$(UTILITIES_DIR)/GenomeFile/ \
           -I$(UTILITIES_DIR)/lineFileUtilities/ \
           -I$(UTILITIES_DIR)/gzstream/ \
           -I$(UTILITIES_DIR)/fileType/ \
           -I$(UTILITIES_DIR)/general/ \
           -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/src \
           -I$(UTILITIES_DIR)/version/ \
           -I$(UTILITIES_DIR)/FileRecordTools/ \
		   -I$(UTILITIES_DIR)/FileRecordTools/FileReaders \
           -I$(UTILITIES_DIR)/FileRecordTools/Records \
  			-I$(UTILITIES_DIR)/KeyListOps/ \
           -I$(UTILITIES_DIR)/general

# ----------------------------------
# define our source and object files
# ----------------------------------
#SOURCES= regressTestMain.cpp RegressTest.cpp RegressTest.h
#OBJECTS= regressTestMain.o RegressTest.o
#BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= regressTest


all: #$(BUILT_OBJECTS)
	@echo "compiling RegressTest.cpp"
	@$(CXX) -c -o $(OBJ_DIR)/RegressTest.o RegressTest.cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)
	@echo "compiling regressTestMain.cpp"
	@$(CXX) -c -o $(OBJ_DIR)/regressTestMain.o regressTestMain.cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)
	
.PHONY: all

#$(BUILT_OBJECTS): $(SOURCES)
#	@echo "  * compiling" $(*F).cpp
#	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)


clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/regressTestMain.o $(OBJ_DIR)/regressTest.o

.PHONY: clean
