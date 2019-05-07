BUILT_OBJECTS += obj/ToolBase.o

obj/ToolBase.o: src/utils/ToolBase/ToolBase.cpp obj/ToolBase.d
	$(CXX_COMPILE)
