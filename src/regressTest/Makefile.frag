BUILT_OBJECTS += obj/regressTestMain.o obj/RegressTest.o

obj/regressTestMain.o: src/regressTest/regressTestMain.cpp obj/regressTestMain.d
	$(CXX_COMPILE)

obj/RegressTest.o: src/regressTest/RegressTest.cpp obj/RegressTest.d
	$(CXX_COMPILE)
