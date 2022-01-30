src/utils/BamTools/include/BamAlignment.mapping.hpp: src/utils/BamTools/mapping/BamAlignment.py src/utils/BamTools/mapping/BamAlignment.map

src/utils/BamTools/include/%.mapping.hpp: src/utils/BamTools/mapping/%.py src/utils/BamTools/mapping/%.map
	$(PYTHON) $^ > $@
