bin_PROGRAMS = bedtools

bedtools_SOURCES = \
	src/bedtools.cpp

bedtools_LDADD = libbedtools.a $(ZLIB_LIBS)

noinst_LIBRARIES = libbedtools.a

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/closestFile \
	-I$(top_srcdir)/src/complementFile \
	-I$(top_srcdir)/src/coverageFile \
	-I$(top_srcdir)/src/fisher \
	-I$(top_srcdir)/src/groupBy \
	-I$(top_srcdir)/src/intersectFile \
	-I$(top_srcdir)/src/jaccard \
	-I$(top_srcdir)/src/mapFile \
	-I$(top_srcdir)/src/mergeFile \
	-I$(top_srcdir)/src/sampleFile \
	-I$(top_srcdir)/src/spacingFile \
	-I$(top_srcdir)/src/subtractFile \
	-I$(top_srcdir)/src/utils/BamTools \
	-I$(top_srcdir)/src/utils/BamTools-Ancillary \
	-I$(top_srcdir)/src/utils/BamTools/include \
	-I$(top_srcdir)/src/utils/BamTools/src \
	-I$(top_srcdir)/src/utils/bedFile \
	-I$(top_srcdir)/src/utils/bedFilePE \
	-I$(top_srcdir)/src/utils/bedGraphFile \
	-I$(top_srcdir)/src/utils/BinTree \
	-I$(top_srcdir)/src/utils/BlockedIntervals \
	-I$(top_srcdir)/src/utils/chromsweep \
	-I$(top_srcdir)/src/utils/Contexts \
	-I$(top_srcdir)/src/utils/driver \
	-I$(top_srcdir)/src/utils/Fasta \
	-I$(top_srcdir)/src/utils/FileRecordTools \
	-I$(top_srcdir)/src/utils/FileRecordTools/FileReaders \
	-I$(top_srcdir)/src/utils/FileRecordTools/Records \
	-I$(top_srcdir)/src/utils/fileType \
	-I$(top_srcdir)/src/utils/general \
	-I$(top_srcdir)/src/utils/GenomeFile \
	-I$(top_srcdir)/src/utils/gzstream \
	-I$(top_srcdir)/src/utils/KeyListOps \
	-I$(top_srcdir)/src/utils/lineFileUtilities \
	-I$(top_srcdir)/src/utils/NewChromsweep \
	-I$(top_srcdir)/src/utils/Point \
	-I$(top_srcdir)/src/utils/RecordOutputMgr \
	-I$(top_srcdir)/src/utils/sequenceUtilities \
	-I$(top_srcdir)/src/utils/stringUtilities \
	-I$(top_srcdir)/src/utils/tabFile \
	-I$(top_srcdir)/src/utils/ToolBase \
	-I$(top_srcdir)/src/utils/VectorOps \
	-I$(top_srcdir)/src/utils/version \
	$(ZLIB_CFLAGS)

libbedtools_a_SOURCES = \
	src/annotateBed/annotateBed.cpp \
	src/annotateBed/annotateBed.h \
	src/annotateBed/annotateMain.cpp \
	src/bamToBed/bamToBed.cpp \
	src/bamToFastq/bamToFastq.cpp \
	src/bamToFastq/bamToFastq.h \
	src/bamToFastq/bamToFastqMain.cpp \
	src/bed12ToBed6/bed12ToBed6.cpp \
	src/bedpeToBam/bedpeToBam.cpp \
	src/bedToBam/bedToBam.cpp \
	src/bedToIgv/bedToIgv.cpp \
	src/closestFile/closestFile.cpp \
	src/closestFile/closestFile.h \
	src/closestFile/closestHelp.cpp \
	src/clusterBed/clusterBed.cpp \
	src/clusterBed/clusterBed.h \
	src/clusterBed/clusterMain.cpp \
	src/complementFile/complementFile.cpp \
	src/complementFile/complementFile.h \
	src/complementFile/complementHelp.cpp \
	src/coverageFile/coverageFile.cpp \
	src/coverageFile/coverageFile.h \
	src/coverageFile/coverageHelp.cpp \
	src/expand/expand.cpp \
	src/fastaFromBed/fastaFromBed.cpp \
	src/fastaFromBed/fastaFromBed.h \
	src/fastaFromBed/fastaFromBedMain.cpp \
	src/fisher/fisher.cpp \
	src/fisher/fisher.h \
	src/fisher/fisherHelp.cpp \
	src/fisher/kfunc.cpp \
	src/fisher/kfunc.h \
	src/flankBed/flankBed.cpp \
	src/flankBed/flankBed.h \
	src/flankBed/flankBedMain.cpp \
	src/genomeCoverageBed/genomeCoverageBed.cpp \
	src/genomeCoverageBed/genomeCoverageBed.h \
	src/genomeCoverageBed/genomeCoverageMain.cpp \
	src/getOverlap/getOverlap.cpp \
	src/groupBy/groupBy.cpp \
	src/groupBy/groupBy.h \
	src/groupBy/groupByHelp.cpp \
	src/intersectFile/intersectFile.cpp \
	src/intersectFile/intersectFile.h \
	src/intersectFile/intersectHelp.cpp \
	src/jaccard/jaccard.cpp \
	src/jaccard/jaccard.h \
	src/jaccard/jaccardHelp.cpp \
	src/linksBed/linksBed.cpp \
	src/linksBed/linksBed.h \
	src/linksBed/linksMain.cpp \
	src/mapFile/mapFile.cpp \
	src/mapFile/mapFile.h \
	src/mapFile/mapHelp.cpp \
	src/maskFastaFromBed/maskFastaFromBed.cpp \
	src/maskFastaFromBed/maskFastaFromBed.h \
	src/maskFastaFromBed/maskFastaFromBedMain.cpp \
	src/mergeFile/mergeFile.cpp \
	src/mergeFile/mergeFile.h \
	src/mergeFile/mergeHelp.cpp \
	src/multiBamCov/multiBamCov.cpp \
	src/multiBamCov/multiBamCov.h \
	src/multiBamCov/multiBamCovMain.cpp \
	src/multiIntersectBed/multiIntersectBed.cpp \
	src/multiIntersectBed/multiIntersectBed.h \
	src/multiIntersectBed/multiIntersectBedMain.cpp \
	src/nekSandbox1/nekSandboxMain.cpp \
	src/nucBed/nucBed.cpp \
	src/nucBed/nucBed.h \
	src/nucBed/nucBedMain.cpp \
	src/pairToBed/pairToBed.cpp \
	src/pairToBed/pairToBed.h \
	src/pairToBed/pairToBedMain.cpp \
	src/pairToPair/pairToPair.cpp \
	src/pairToPair/pairToPair.h \
	src/pairToPair/pairToPairMain.cpp \
	src/randomBed/randomBed.cpp \
	src/randomBed/randomBed.h \
	src/randomBed/randomBedMain.cpp \
	src/regressTest/RegressTest.cpp \
	src/regressTest/RegressTest.h \
	src/regressTest/regressTestMain.cpp \
	src/reldist/reldist.cpp \
	src/reldist/reldist.h \
	src/reldist/reldistMain.cpp \
	src/sampleFile/sampleFile.cpp \
	src/sampleFile/sampleFile.h \
	src/sampleFile/sampleHelp.cpp \
	src/shiftBed/shiftBed.cpp \
	src/shiftBed/shiftBed.h \
	src/shiftBed/shiftBedMain.cpp \
	src/shuffleBed/shuffleBed.cpp \
	src/shuffleBed/shuffleBed.h \
	src/shuffleBed/shuffleBedMain.cpp \
	src/slopBed/slopBed.cpp \
	src/slopBed/slopBed.h \
	src/slopBed/slopBedMain.cpp \
	src/sortBed/sortBed.cpp \
	src/sortBed/sortBed.h \
	src/sortBed/sortMain.cpp \
	src/spacingFile/spacingFile.cpp \
	src/spacingFile/spacingFile.h \
	src/spacingFile/spacingHelp.cpp \
	src/split/splitBed.cpp \
	src/split/splitBed.h \
	src/split/splitBedMain.cpp \
	src/subtractFile/subtractFile.cpp \
	src/subtractFile/subtractFile.h \
	src/subtractFile/subtractHelp.cpp \
	src/tagBam/tagBam.cpp \
	src/tagBam/tagBam.h \
	src/tagBam/tagBamMain.cpp \
	src/unionBedGraphs/unionBedGraphs.cpp \
	src/unionBedGraphs/unionBedGraphs.h \
	src/unionBedGraphs/unionBedGraphsMain.cpp \
	src/utils/BamTools-Ancillary/BamAncillary.cpp \
	src/utils/BamTools-Ancillary/BamAncillary.h \
	src/utils/bedFile/bedFile.cpp \
	src/utils/bedFile/bedFile.h \
	src/utils/bedFilePE/bedFilePE.cpp \
	src/utils/bedFilePE/bedFilePE.h \
	src/utils/bedGraphFile/bedGraphFile.cpp \
	src/utils/bedGraphFile/bedGraphFile.h \
	src/utils/BinTree/BinTree.cpp \
	src/utils/BinTree/BinTree.h \
	src/utils/BlockedIntervals/BlockedIntervals.cpp \
	src/utils/BlockedIntervals/BlockedIntervals.h \
	src/utils/chromsweep/chromsweep.cpp \
	src/utils/chromsweep/chromsweep.h \
	src/utils/Contexts/ContextBase.cpp \
	src/utils/Contexts/ContextBase.h \
	src/utils/Contexts/ContextClosest.cpp \
	src/utils/Contexts/ContextClosest.h \
	src/utils/Contexts/ContextComplement.cpp \
	src/utils/Contexts/ContextComplement.h \
	src/utils/Contexts/ContextCoverage.cpp \
	src/utils/Contexts/ContextCoverage.h \
	src/utils/Contexts/ContextFisher.cpp \
	src/utils/Contexts/ContextFisher.h \
	src/utils/Contexts/ContextGroupBy.cpp \
	src/utils/Contexts/ContextGroupBy.h \
	src/utils/Contexts/ContextIntersect.cpp \
	src/utils/Contexts/ContextIntersect.h \
	src/utils/Contexts/ContextJaccard.cpp \
	src/utils/Contexts/ContextJaccard.h \
	src/utils/Contexts/ContextMap.cpp \
	src/utils/Contexts/ContextMap.h \
	src/utils/Contexts/ContextMerge.cpp \
	src/utils/Contexts/ContextMerge.h \
	src/utils/Contexts/ContextSample.cpp \
	src/utils/Contexts/ContextSample.h \
	src/utils/Contexts/ContextSpacing.cpp \
	src/utils/Contexts/ContextSpacing.h \
	src/utils/Contexts/ContextSubtract.cpp \
	src/utils/Contexts/ContextSubtract.h \
	src/utils/driver/BedtoolsDriver.cpp \
	src/utils/driver/BedtoolsDriver.h \
	src/utils/Fasta/Fasta.cpp \
	src/utils/Fasta/Fasta.h \
	src/utils/Fasta/LargeFileSupport.h \
	src/utils/Fasta/split.cpp \
	src/utils/Fasta/split.h \
	src/utils/FileRecordTools/FileReaders/BamFileReader.cpp \
	src/utils/FileRecordTools/FileReaders/BamFileReader.h \
	src/utils/FileRecordTools/FileReaders/BufferedStreamMgr.cpp \
	src/utils/FileRecordTools/FileReaders/BufferedStreamMgr.h \
	src/utils/FileRecordTools/FileReaders/FileReader.cpp \
	src/utils/FileRecordTools/FileReaders/FileReader.h \
	src/utils/FileRecordTools/FileReaders/InputStreamMgr.cpp \
	src/utils/FileRecordTools/FileReaders/InputStreamMgr.h \
	src/utils/FileRecordTools/FileReaders/SingleLineDelimTextFileReader.cpp \
	src/utils/FileRecordTools/FileReaders/SingleLineDelimTextFileReader.h \
	src/utils/FileRecordTools/FileReaders/SingleLineDelimTransferBuffer.cpp \
	src/utils/FileRecordTools/FileReaders/SingleLineDelimTransferBuffer.h \
	src/utils/FileRecordTools/FileRecordMergeMgr.cpp \
	src/utils/FileRecordTools/FileRecordMergeMgr.h \
	src/utils/FileRecordTools/FileRecordMgr.cpp \
	src/utils/FileRecordTools/FileRecordMgr.h \
	src/utils/FileRecordTools/Records/BamRecord.cpp \
	src/utils/FileRecordTools/Records/BamRecord.h \
	src/utils/FileRecordTools/Records/Bed12Interval.cpp \
	src/utils/FileRecordTools/Records/Bed12Interval.h \
	src/utils/FileRecordTools/Records/Bed3Interval.cpp \
	src/utils/FileRecordTools/Records/Bed3Interval.h \
	src/utils/FileRecordTools/Records/Bed4Interval.cpp \
	src/utils/FileRecordTools/Records/Bed4Interval.h \
	src/utils/FileRecordTools/Records/Bed5Interval.cpp \
	src/utils/FileRecordTools/Records/Bed5Interval.h \
	src/utils/FileRecordTools/Records/Bed6Interval.cpp \
	src/utils/FileRecordTools/Records/Bed6Interval.h \
	src/utils/FileRecordTools/Records/BedGraphInterval.cpp \
	src/utils/FileRecordTools/Records/BedGraphInterval.h \
	src/utils/FileRecordTools/Records/BedPlusInterval.cpp \
	src/utils/FileRecordTools/Records/BedPlusInterval.h \
	src/utils/FileRecordTools/Records/BlockMgr.cpp \
	src/utils/FileRecordTools/Records/BlockMgr.h \
	src/utils/FileRecordTools/Records/EmptyRecord.cpp \
	src/utils/FileRecordTools/Records/EmptyRecord.h \
	src/utils/FileRecordTools/Records/GffPlusRecord.cpp \
	src/utils/FileRecordTools/Records/GffPlusRecord.h \
	src/utils/FileRecordTools/Records/GffRecord.cpp \
	src/utils/FileRecordTools/Records/GffRecord.h \
	src/utils/FileRecordTools/Records/NoPosPlusRecord.cpp \
	src/utils/FileRecordTools/Records/NoPosPlusRecord.h \
	src/utils/FileRecordTools/Records/PlusFields.cpp \
	src/utils/FileRecordTools/Records/PlusFields.h \
	src/utils/FileRecordTools/Records/Record.cpp \
	src/utils/FileRecordTools/Records/Record.h \
	src/utils/FileRecordTools/Records/RecordKeyList.cpp \
	src/utils/FileRecordTools/Records/RecordKeyList.h \
	src/utils/FileRecordTools/Records/RecordKeyVector.cpp \
	src/utils/FileRecordTools/Records/RecordKeyVector.h \
	src/utils/FileRecordTools/Records/RecordList.cpp \
	src/utils/FileRecordTools/Records/RecordList.h \
	src/utils/FileRecordTools/Records/RecordMgr.cpp \
	src/utils/FileRecordTools/Records/RecordMgr.h \
	src/utils/FileRecordTools/Records/StrandQueue.cpp \
	src/utils/FileRecordTools/Records/StrandQueue.h \
	src/utils/FileRecordTools/Records/VcfRecord.cpp \
	src/utils/FileRecordTools/Records/VcfRecord.h \
	src/utils/fileType/FileRecordTypeChecker.cpp \
	src/utils/fileType/FileRecordTypeChecker.h \
	src/utils/fileType/fileType.cpp \
	src/utils/fileType/fileType.h \
	src/utils/general/BedtoolsTypes.h \
	src/utils/general/BTlist.h \
	src/utils/general/CommonHelp.cpp \
	src/utils/general/CommonHelp.h \
	src/utils/general/CompressionTools.cpp \
	src/utils/general/CompressionTools.h \
	src/utils/general/ErrorMsg.cpp \
	src/utils/general/ErrorMsg.h \
	src/utils/general/FreeList.h \
	src/utils/general/InflateStreamBuf.h \
	src/utils/general/ParseTools.cpp \
	src/utils/general/ParseTools.h \
	src/utils/general/PushBackStreamBuf.cpp \
	src/utils/general/PushBackStreamBuf.h \
	src/utils/general/Tokenizer.cpp \
	src/utils/general/Tokenizer.h \
	src/utils/GenomeFile/GenomeFile.cpp \
	src/utils/GenomeFile/GenomeFile.h \
	src/utils/GenomeFile/NewGenomeFile.cpp \
	src/utils/GenomeFile/NewGenomeFile.h \
	src/utils/gzstream/gzstream.C \
	src/utils/gzstream/gzstream.h \
	src/utils/KeyListOps/KeyListOps.cpp \
	src/utils/KeyListOps/KeyListOps.h \
	src/utils/KeyListOps/KeyListOpsMethods.cpp \
	src/utils/KeyListOps/KeyListOpsMethods.h \
	src/utils/lineFileUtilities/lineFileUtilities.h \
	src/utils/NewChromsweep/CloseSweep.cpp \
	src/utils/NewChromsweep/CloseSweep.h \
	src/utils/NewChromsweep/NewChromsweep.cpp \
	src/utils/NewChromsweep/NewChromsweep.h \
	src/utils/Point/Point.h \
	src/utils/RecordOutputMgr/RecordOutputMgr.cpp \
	src/utils/RecordOutputMgr/RecordOutputMgr.h \
	src/utils/sequenceUtilities/sequenceUtils.cpp \
	src/utils/sequenceUtilities/sequenceUtils.h \
	src/utils/stringUtilities/stringUtilities.h \
	src/utils/tabFile/tabFile.cpp \
	src/utils/tabFile/tabFile.h \
	src/utils/ToolBase/ToolBase.cpp \
	src/utils/ToolBase/ToolBase.h \
	src/utils/VectorOps/VectorOps.cpp \
	src/utils/VectorOps/VectorOps.h \
	src/windowBed/windowBed.cpp \
	src/windowBed/windowBed.h \
	src/windowBed/windowMain.cpp \
	src/windowMaker/windowMaker.cpp \
	src/windowMaker/windowMaker.h \
	src/windowMaker/windowMakerMain.cpp

if BUNDLED_BAMTOOLS
libbedtools_a_SOURCES += \
	src/utils/BamTools/src/api/algorithms/Sort.h \
	src/utils/BamTools/src/api/api_global.h \
	src/utils/BamTools/src/api/BamAlignment.cpp \
	src/utils/BamTools/src/api/BamAlignment.h \
	src/utils/BamTools/src/api/BamAux.h \
	src/utils/BamTools/src/api/BamConstants.h \
	src/utils/BamTools/src/api/BamIndex.h \
	src/utils/BamTools/src/api/BamMultiReader.cpp \
	src/utils/BamTools/src/api/BamMultiReader.h \
	src/utils/BamTools/src/api/BamReader.cpp \
	src/utils/BamTools/src/api/BamReader.h \
	src/utils/BamTools/src/api/BamWriter.cpp \
	src/utils/BamTools/src/api/BamWriter.h \
	src/utils/BamTools/src/api/IBamIODevice.h \
	src/utils/BamTools/src/api/internal/bam/BamHeader_p.cpp \
	src/utils/BamTools/src/api/internal/bam/BamHeader_p.h \
	src/utils/BamTools/src/api/internal/bam/BamMultiMerger_p.h \
	src/utils/BamTools/src/api/internal/bam/BamMultiReader_p.cpp \
	src/utils/BamTools/src/api/internal/bam/BamMultiReader_p.h \
	src/utils/BamTools/src/api/internal/bam/BamRandomAccessController_p.cpp \
	src/utils/BamTools/src/api/internal/bam/BamRandomAccessController_p.h \
	src/utils/BamTools/src/api/internal/bam/BamReader_p.cpp \
	src/utils/BamTools/src/api/internal/bam/BamReader_p.h \
	src/utils/BamTools/src/api/internal/bam/BamWriter_p.cpp \
	src/utils/BamTools/src/api/internal/bam/BamWriter_p.h \
	src/utils/BamTools/src/api/internal/index/BamIndexFactory_p.cpp \
	src/utils/BamTools/src/api/internal/index/BamIndexFactory_p.h \
	src/utils/BamTools/src/api/internal/index/BamStandardIndex_p.cpp \
	src/utils/BamTools/src/api/internal/index/BamStandardIndex_p.h \
	src/utils/BamTools/src/api/internal/index/BamToolsIndex_p.cpp \
	src/utils/BamTools/src/api/internal/index/BamToolsIndex_p.h \
	src/utils/BamTools/src/api/internal/io/BamDeviceFactory_p.cpp \
	src/utils/BamTools/src/api/internal/io/BamDeviceFactory_p.h \
	src/utils/BamTools/src/api/internal/io/BamFile_p.cpp \
	src/utils/BamTools/src/api/internal/io/BamFile_p.h \
	src/utils/BamTools/src/api/internal/io/BamFtp_p.cpp \
	src/utils/BamTools/src/api/internal/io/BamFtp_p.h \
	src/utils/BamTools/src/api/internal/io/BamHttp_p.cpp \
	src/utils/BamTools/src/api/internal/io/BamHttp_p.h \
	src/utils/BamTools/src/api/internal/io/BamPipe_p.cpp \
	src/utils/BamTools/src/api/internal/io/BamPipe_p.h \
	src/utils/BamTools/src/api/internal/io/BgzfStream_p.cpp \
	src/utils/BamTools/src/api/internal/io/BgzfStream_p.h \
	src/utils/BamTools/src/api/internal/io/ByteArray_p.cpp \
	src/utils/BamTools/src/api/internal/io/ByteArray_p.h \
	src/utils/BamTools/src/api/internal/io/HostAddress_p.cpp \
	src/utils/BamTools/src/api/internal/io/HostAddress_p.h \
	src/utils/BamTools/src/api/internal/io/HostInfo_p.cpp \
	src/utils/BamTools/src/api/internal/io/HostInfo_p.h \
	src/utils/BamTools/src/api/internal/io/HttpHeader_p.cpp \
	src/utils/BamTools/src/api/internal/io/HttpHeader_p.h \
	src/utils/BamTools/src/api/internal/io/ILocalIODevice_p.cpp \
	src/utils/BamTools/src/api/internal/io/ILocalIODevice_p.h \
	src/utils/BamTools/src/api/internal/io/NetUnix_p.h \
	src/utils/BamTools/src/api/internal/io/RollingBuffer_p.cpp \
	src/utils/BamTools/src/api/internal/io/RollingBuffer_p.h \
	src/utils/BamTools/src/api/internal/io/StdStreamDevice_p.h \
	src/utils/BamTools/src/api/internal/io/TcpSocketEngine_p.cpp \
	src/utils/BamTools/src/api/internal/io/TcpSocketEngine_p.h \
	src/utils/BamTools/src/api/internal/io/TcpSocketEngine_unix_p.cpp \
	src/utils/BamTools/src/api/internal/io/TcpSocket_p.cpp \
	src/utils/BamTools/src/api/internal/io/TcpSocket_p.h \
	src/utils/BamTools/src/api/internal/sam/SamFormatParser_p.cpp \
	src/utils/BamTools/src/api/internal/sam/SamFormatParser_p.h \
	src/utils/BamTools/src/api/internal/sam/SamFormatPrinter_p.cpp \
	src/utils/BamTools/src/api/internal/sam/SamFormatPrinter_p.h \
	src/utils/BamTools/src/api/internal/sam/SamHeaderValidator_p.cpp \
	src/utils/BamTools/src/api/internal/sam/SamHeaderValidator_p.h \
	src/utils/BamTools/src/api/internal/sam/SamHeaderVersion_p.h \
	src/utils/BamTools/src/api/internal/utils/BamException_p.cpp \
	src/utils/BamTools/src/api/internal/utils/BamException_p.h \
	src/utils/BamTools/src/api/SamConstants.h \
	src/utils/BamTools/src/api/SamHeader.cpp \
	src/utils/BamTools/src/api/SamHeader.h \
	src/utils/BamTools/src/api/SamProgramChain.cpp \
	src/utils/BamTools/src/api/SamProgramChain.h \
	src/utils/BamTools/src/api/SamProgram.cpp \
	src/utils/BamTools/src/api/SamProgram.h \
	src/utils/BamTools/src/api/SamReadGroup.cpp \
	src/utils/BamTools/src/api/SamReadGroup.h \
	src/utils/BamTools/src/api/SamReadGroupDictionary.cpp \
	src/utils/BamTools/src/api/SamReadGroupDictionary.h \
	src/utils/BamTools/src/api/SamSequence.cpp \
	src/utils/BamTools/src/api/SamSequence.h \
	src/utils/BamTools/src/api/SamSequenceDictionary.cpp \
	src/utils/BamTools/src/api/SamSequenceDictionary.h \
	src/utils/BamTools/src/shared/bamtools_global.h
else
AM_CPPFLAGS += $(BAMTOOLS_CFLAGS)
bedtools_LDADD += $(BAMTOOLS_LIBS)
endif

if ENABLE_DOCUMENTATION
SUBDIRS = docs
endif

EXTRA_DIST = test scripts/makeBashScripts.py git-version-gen $(top_srcdir)/.version

BUILT_SOURCES = $(top_srcdir)/.version
$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@
dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version

dist_doc_DATA = README.md RELEASE_HISTORY

genomesdir = $(datadir)/bedtools/genomes
dist_genomes_DATA = \
	genomes/mouse.mm8.genome \
	genomes/mouse.mm9.genome \
	genomes/mouse.mm10.genome \
	genomes/human.hg18.genome \
	genomes/human.hg19.genome \
	genomes/human.hg38.genome

# we copy the original test dir 'test'
# to a folder called 'tests', in order
# to circumvent distcheck's write checks.
# After running the testsuite, we can
# safely delete the folder again
check: all
	@echo "Running test suite"
	cp -R $(top_srcdir)/test tests
	chmod -R u+wx tests
	cd tests && sh test.sh
	rm -rf tests
	@echo "Done running test suite"

install-exec-hook:
	@echo "Creating compatibility scripts for old CLI"
	cd $(DESTDIR)$(bindir) && \
		python $(abs_top_srcdir)/scripts/makeBashScripts.py || exit 1
	chmod +x $(DESTDIR)$(bindir)/*
	@echo "Done creating compatibility scripts for old CLI"

# remove compatibility links
# needed in order for 'make distcheck' to succeed
uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/annotateBed
	rm -f $(DESTDIR)$(bindir)/bamToBed
	rm -f $(DESTDIR)$(bindir)/bamToFastq
	rm -f $(DESTDIR)$(bindir)/bed12ToBed6
	rm -f $(DESTDIR)$(bindir)/bedpeToBam
	rm -f $(DESTDIR)$(bindir)/bedToBam
	rm -f $(DESTDIR)$(bindir)/bedToIgv
	rm -f $(DESTDIR)$(bindir)/closestBed
	rm -f $(DESTDIR)$(bindir)/clusterBed
	rm -f $(DESTDIR)$(bindir)/complementBed
	rm -f $(DESTDIR)$(bindir)/coverageBed
	rm -f $(DESTDIR)$(bindir)/expandCols
	rm -f $(DESTDIR)$(bindir)/fastaFromBed
	rm -f $(DESTDIR)$(bindir)/flankBed
	rm -f $(DESTDIR)$(bindir)/genomeCoverageBed
	rm -f $(DESTDIR)$(bindir)/getOverlap
	rm -f $(DESTDIR)$(bindir)/groupBy
	rm -f $(DESTDIR)$(bindir)/intersectBed
	rm -f $(DESTDIR)$(bindir)/linksBed
	rm -f $(DESTDIR)$(bindir)/mapBed
	rm -f $(DESTDIR)$(bindir)/maskFastaFromBed
	rm -f $(DESTDIR)$(bindir)/mergeBed
	rm -f $(DESTDIR)$(bindir)/multiBamCov
	rm -f $(DESTDIR)$(bindir)/multiIntersectBed
	rm -f $(DESTDIR)$(bindir)/nucBed
	rm -f $(DESTDIR)$(bindir)/pairToBed
	rm -f $(DESTDIR)$(bindir)/pairToPair
	rm -f $(DESTDIR)$(bindir)/randomBed
	rm -f $(DESTDIR)$(bindir)/shiftBed
	rm -f $(DESTDIR)$(bindir)/shuffleBed
	rm -f $(DESTDIR)$(bindir)/slopBed
	rm -f $(DESTDIR)$(bindir)/sortBed
	rm -f $(DESTDIR)$(bindir)/subtractBed
	rm -f $(DESTDIR)$(bindir)/tagBam
	rm -f $(DESTDIR)$(bindir)/unionBedGraphs
	rm -f $(DESTDIR)$(bindir)/windowBed
	rm -f $(DESTDIR)$(bindir)/windowMaker
