all: include/BamAlignment.mapping.hpp
include/%.mapping.hpp: mapping/%.py mapping/%.map
	python $^ > $@
