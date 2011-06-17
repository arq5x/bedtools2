# -------------------
# define our includes
# -------------------
OBJ_DIR = ../../../obj/
INCLUDES = -Isrc/ -Iinclude/

# ----------------------------------
# define our source and object files
# ----------------------------------

SOURCES= src/api/BamAlignment.cpp \
         src/api/BamMultiReader.cpp \
         src/api/BamReader.cpp \
         src/api/BamWriter.cpp \
         src/api/SamHeader.cpp \
         src/api/SamProgram.cpp \
         src/api/SamProgramChain.cpp \
         src/api/SamReadGroup.cpp \
         src/api/SamReadGroupDictionary.cpp \
         src/api/SamSequence.cpp \
         src/api/SamSequenceDictionary.cpp \
		 src/api/internal/BamHeader_p.cpp \
		 src/api/internal/BamIndexFactory_p.cpp \
		 src/api/internal/BamMultiReader_p.cpp \
		 src/api/internal/BamRandomAccessController_p.cpp \
		 src/api/internal/BamReader_p.cpp \
		 src/api/internal/BamStandardIndex_p.cpp \
		 src/api/internal/BamToolsIndex_p.cpp \
		 src/api/internal/BamWriter_p.cpp \
		 src/api/internal/BgzfStream_p.cpp \
		 src/api/internal/SamFormatParser_p.cpp \
		 src/api/internal/SamFormatPrinter_p.cpp \
		 src/api/internal/SamHeaderValidator_p.cpp

# ----------------------------------
# define our source and object files
# ----------------------------------
OBJECTS= $(SOURCES:.cpp=.o)
LIBRARY=libbamtools.a

all: $(LIBRARY)
    
.PHONY: all

$(LIBRARY): $(OBJECTS)
	[ -d lib ] || mkdir -p lib
	@echo "  * linking $(LIBRARY)"
	ar cr lib/$@  $^
		
$(OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c $(*D)/$(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o $(*D)/$(*F).o
	
clean:
	@echo "Cleaning up."
	@rm -f lib/*
	@rm -f src/api/*.o
	@rm -f src/api/internal/*.o

.PHONY: clean