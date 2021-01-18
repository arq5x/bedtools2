BUILT_OBJECTS += obj/ParseTools.o obj/PushBackStreamBuf.o obj/CompressionTools.o obj/Tokenizer.o obj/CommonHelp.o obj/ErrorMsg.o obj/Random.o

obj/ParseTools.o: src/utils/general/ParseTools.cpp obj/ParseTools.d
	$(CXX_COMPILE)

obj/PushBackStreamBuf.o: src/utils/general/PushBackStreamBuf.cpp obj/PushBackStreamBuf.d
	$(CXX_COMPILE)

obj/CompressionTools.o: src/utils/general/CompressionTools.cpp obj/CompressionTools.d
	$(CXX_COMPILE)

obj/Tokenizer.o: src/utils/general/Tokenizer.cpp obj/Tokenizer.d
	$(CXX_COMPILE)

obj/CommonHelp.o: src/utils/general/CommonHelp.cpp obj/CommonHelp.d
	$(CXX_COMPILE)

obj/ErrorMsg.o: src/utils/general/ErrorMsg.cpp obj/ErrorMsg.d
	$(CXX_COMPILE)

obj/Random.o: src/utils/general/Random.cpp obj/Random.d
	$(CXX_COMPILE)
