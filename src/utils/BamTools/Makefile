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
         src/api/internal/bam/BamHeader_p.cpp \
         src/api/internal/bam/BamMultiReader_p.cpp \
         src/api/internal/bam/BamRandomAccessController_p.cpp \
         src/api/internal/bam/BamReader_p.cpp \
         src/api/internal/bam/BamWriter_p.cpp \
         src/api/internal/index/BamIndexFactory_p.cpp \
         src/api/internal/index/BamStandardIndex_p.cpp \
         src/api/internal/index/BamToolsIndex_p.cpp \
         src/api/internal/io/BamDeviceFactory_p.cpp \
         src/api/internal/io/BamFile_p.cpp \
         src/api/internal/io/BamFtp_p.cpp \
         src/api/internal/io/BamHttp_p.cpp \
         src/api/internal/io/BamPipe_p.cpp \
         src/api/internal/io/BgzfStream_p.cpp \
         src/api/internal/io/ByteArray_p.cpp \
         src/api/internal/io/HostAddress_p.cpp \
         src/api/internal/io/HostInfo_p.cpp \
         src/api/internal/io/HttpHeader_p.cpp \
         src/api/internal/io/ILocalIODevice_p.cpp \
         src/api/internal/io/RollingBuffer_p.cpp \
         src/api/internal/io/TcpSocketEngine_p.cpp \
         src/api/internal/io/TcpSocketEngine_unix_p.cpp \
         src/api/internal/io/TcpSocket_p.cpp \
         src/api/internal/sam/SamFormatParser_p.cpp \
         src/api/internal/sam/SamFormatPrinter_p.cpp \
         src/api/internal/sam/SamHeaderValidator_p.cpp \
         src/api/internal/utils/BamException_p.cpp

OBJECTS= $(SOURCES:.cpp=.o)
LIBRARY=lib/libbamtools.a

all: $(LIBRARY)
    
.PHONY: all

$(LIBRARY): $(OBJECTS)
	@echo " * linking $(LIBRARY)"
	ar cr $@ $^
	-ranlib $@

$(OBJECTS): $(SOURCES)
	@echo " * compiling" $(*F).cpp
	@$(CXX) -c $(*D)/$(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o $(*D)/$(*F).o

.PHONY: api
api:
	@echo " * Creating BamTools API"
	@mkdir -p lib
	@mkdir -p include
	@mkdir -p include/api
	@mkdir -p include/shared
	@cp src/api/*.h include/api
	@cp src/shared/*.h include/shared

.PHONY: clean_api
clean_api:
	@echo " * Cleaning-up BamTools API"
	@rm -rf ./lib ./src/api/*.o ./src/api/internal/*/*.o ./include
