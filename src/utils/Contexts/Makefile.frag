BUILT_OBJECTS += obj/ContextBase.o obj/ContextIntersect.o obj/ContextFisher.o obj/ContextMap.o obj/ContextSample.o obj/ContextSpacing.o obj/ContextMerge.o obj/ContextJaccard.o obj/ContextClosest.o obj/ContextSubtract.o obj/ContextCoverage.o obj/ContextComplement.o obj/ContextGroupBy.o obj/ContextSummary.o

obj/ContextBase.o: src/utils/Contexts/ContextBase.cpp obj/ContextBase.d
	$(CXX_COMPILE)

obj/ContextIntersect.o: src/utils/Contexts/ContextIntersect.cpp obj/ContextIntersect.d
	$(CXX_COMPILE)

obj/ContextFisher.o: src/utils/Contexts/ContextFisher.cpp obj/ContextFisher.d
	$(CXX_COMPILE)

obj/ContextMap.o: src/utils/Contexts/ContextMap.cpp obj/ContextMap.d
	$(CXX_COMPILE)

obj/ContextSample.o: src/utils/Contexts/ContextSample.cpp obj/ContextSample.d
	$(CXX_COMPILE)

obj/ContextSpacing.o: src/utils/Contexts/ContextSpacing.cpp obj/ContextSpacing.d
	$(CXX_COMPILE)

obj/ContextMerge.o: src/utils/Contexts/ContextMerge.cpp obj/ContextMerge.d
	$(CXX_COMPILE)

obj/ContextJaccard.o: src/utils/Contexts/ContextJaccard.cpp obj/ContextJaccard.d
	$(CXX_COMPILE)

obj/ContextClosest.o: src/utils/Contexts/ContextClosest.cpp obj/ContextClosest.d
	$(CXX_COMPILE)

obj/ContextSubtract.o: src/utils/Contexts/ContextSubtract.cpp obj/ContextSubtract.d
	$(CXX_COMPILE)

obj/ContextCoverage.o: src/utils/Contexts/ContextCoverage.cpp obj/ContextCoverage.d
	$(CXX_COMPILE)

obj/ContextComplement.o: src/utils/Contexts/ContextComplement.cpp obj/ContextComplement.d
	$(CXX_COMPILE)

obj/ContextGroupBy.o: src/utils/Contexts/ContextGroupBy.cpp obj/ContextGroupBy.d
	$(CXX_COMPILE)

obj/ContextSummary.o: src/utils/Contexts/ContextSummary.cpp obj/ContextSummary.d
	$(CXX_COMPILE)
