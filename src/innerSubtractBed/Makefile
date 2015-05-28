UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= innerSubtractBed.pl
BUILT_OBJECTS= $(patsubst %.pl,$(BIN_DIR)/%,$(SOURCES))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(@F)
	@cp $(@F).pl $@
	
clean:
	@echo "Cleaning up."

.PHONY: clean