#include "RecordMgr.h"

#include "FreeList.h"

RecordMgr::RecordMgr(FileRecordTypeChecker::RECORD_TYPE recType, int blockSize)
: _recordType(recType),
  _freeList(NULL),
  _freeListBlockSize(blockSize)
{
	switch(_recordType) {
		case FileRecordTypeChecker::EMPTY_RECORD_TYPE:
		{
			_freeList = new FreeList<EmptyRecord>(0);
			break;
		}


		case FileRecordTypeChecker::BED3_RECORD_TYPE:
		{
			_freeList = new FreeList<Bed3Interval>(_freeListBlockSize);
			break;
		}

		case FileRecordTypeChecker::BED4_RECORD_TYPE:
		{
			_freeList = new FreeList<Bed4Interval>(_freeListBlockSize);
			break;
		}

		case FileRecordTypeChecker::BED5_RECORD_TYPE:
		{
			_freeList = new FreeList<Bed5Interval>(_freeListBlockSize);
			break;
		}

		case FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE:
		{
			_freeList = new FreeList<BedGraphInterval>(_freeListBlockSize);
			break;
		}

		case FileRecordTypeChecker::BED6_RECORD_TYPE:
		{
			_freeList = new FreeList<Bed6Interval>(_freeListBlockSize);
			break;
		}
		case FileRecordTypeChecker::BED_PLUS_RECORD_TYPE:
		case FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE:
		{
			_freeList = new FreeList<BedPlusInterval>(_freeListBlockSize);
			break;
		}
		case FileRecordTypeChecker::BED12_RECORD_TYPE:
		{
			_freeList = new FreeList<Bed12Interval>(_freeListBlockSize);
			break;
		}
		case FileRecordTypeChecker::BAM_RECORD_TYPE:
		{
			_freeList = new FreeList<BamRecord>(_freeListBlockSize);
			break;
		}
		case FileRecordTypeChecker::VCF_RECORD_TYPE:
		{
			_freeList = new FreeList<VcfRecord>(_freeListBlockSize);
			break;
		}
		case FileRecordTypeChecker::GFF_RECORD_TYPE:
		{
			_freeList = new FreeList<GffRecord>(_freeListBlockSize);
			break;
		}
		case FileRecordTypeChecker::GFF_PLUS_RECORD_TYPE:
		{
			_freeList = new FreeList<GffPlusRecord>(_freeListBlockSize);
			break;
		}
		case FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE:
		{
			_freeList = new FreeList<NoPosPlusRecord>(_freeListBlockSize);
			break;
		}


		default:
			break;
	}
}

RecordMgr::~RecordMgr()
{
	switch(_recordType) {
		case FileRecordTypeChecker::EMPTY_RECORD_TYPE:
		{
			delete (FreeList<EmptyRecord> *)_freeList;
			break;
		}
		case FileRecordTypeChecker::BED3_RECORD_TYPE:
		{
			delete (FreeList<Bed3Interval> *)_freeList;
			break;
		}

		case FileRecordTypeChecker::BED4_RECORD_TYPE:
		{
			delete (FreeList<Bed4Interval> *)_freeList;
			break;
		}

		case FileRecordTypeChecker::BED5_RECORD_TYPE:
		{
			delete (FreeList<Bed5Interval> *)_freeList;
			break;
		}

		case FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE:
		{
			delete (FreeList<BedGraphInterval> *)_freeList;
			break;
		}

		case FileRecordTypeChecker::BED6_RECORD_TYPE:
		{
			delete (FreeList<Bed6Interval> *)_freeList;
			break;
		}
		case FileRecordTypeChecker::BED_PLUS_RECORD_TYPE:
		case FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE:
		{
			delete (FreeList<BedPlusInterval> *)_freeList;
			break;
		}
		case FileRecordTypeChecker::BED12_RECORD_TYPE:
		{
			delete (FreeList<Bed12Interval> *)_freeList;
			break;
		}
		case FileRecordTypeChecker::BAM_RECORD_TYPE:
		{
			delete (FreeList<BamRecord> *)_freeList;
			break;
		}
		case FileRecordTypeChecker::VCF_RECORD_TYPE:
		{
			delete (FreeList<VcfRecord> *)_freeList;
			break;
		}
		case FileRecordTypeChecker::GFF_RECORD_TYPE:
		{
			delete (FreeList<GffRecord> *)_freeList;
			break;
		}
		case FileRecordTypeChecker::GFF_PLUS_RECORD_TYPE:
		{
			delete (FreeList<GffPlusRecord> *)_freeList;
			break;
		}
		case FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE:
		{
			delete (FreeList<NoPosPlusRecord> *)_freeList;
			break;
		}


		default:
			break;
	}
}

Record *RecordMgr::allocateRecord()
{
	Record *record = NULL;
	switch(_recordType) {
		case FileRecordTypeChecker::EMPTY_RECORD_TYPE:
		{
			EmptyRecord *er = ((FreeList<EmptyRecord> *)_freeList)->newObj();
			record = er;
			break;
		}
		case FileRecordTypeChecker::BED3_RECORD_TYPE:
		{
			Bed3Interval *b3i = ((FreeList<Bed3Interval> *)_freeList)->newObj();
			record = b3i;
			break;
		}

		case FileRecordTypeChecker::BED4_RECORD_TYPE:
		{
			Bed4Interval *b4i = ((FreeList<Bed4Interval> *)_freeList)->newObj();
			record = b4i;
			break;
		}

		case FileRecordTypeChecker::BED5_RECORD_TYPE:
		{
			Bed5Interval *b5i = ((FreeList<Bed5Interval> *)_freeList)->newObj();
			record = b5i;
			break;
		}

		case FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE:
		{
			BedGraphInterval *bgi = ((FreeList<BedGraphInterval> *)_freeList)->newObj();
			record = bgi;
			break;
		}

		case FileRecordTypeChecker::BED6_RECORD_TYPE:
		{
			Bed6Interval *b6i = ((FreeList<Bed6Interval> *)_freeList)->newObj();
			record = b6i;
			break;
		}
		case FileRecordTypeChecker::BED_PLUS_RECORD_TYPE:
		case FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE:
		{
			BedPlusInterval *bPi = ((FreeList<BedPlusInterval> *)_freeList)->newObj();
			if (_recordType == FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE) {
				bPi->setNumFixedFields(6);
			}
			record = bPi;
			break;
		}
		case FileRecordTypeChecker::BED12_RECORD_TYPE:
		{
			Bed12Interval *b12i = ((FreeList<Bed12Interval> *)_freeList)->newObj();
			record = b12i;
			break;
		}
		case FileRecordTypeChecker::BAM_RECORD_TYPE:
		{
			BamRecord *bamRec = ((FreeList<BamRecord> *)_freeList)->newObj();
			record = bamRec;
			break;
		}

		case FileRecordTypeChecker::VCF_RECORD_TYPE:
		{
			VcfRecord *vcfRec = ((FreeList<VcfRecord> *)_freeList)->newObj();
			record = vcfRec;
			break;
		}

		case FileRecordTypeChecker::GFF_RECORD_TYPE:
		{
			GffRecord *gfr = ((FreeList<GffRecord> *)_freeList)->newObj();
			record = gfr;
			break;
		}
		case FileRecordTypeChecker::GFF_PLUS_RECORD_TYPE:
		{
			GffPlusRecord *gfpr = ((FreeList<GffPlusRecord> *)_freeList)->newObj();
			record = gfpr;
			break;
		}
		case FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE:
		{
			NoPosPlusRecord *nppr = ((FreeList<NoPosPlusRecord> *)_freeList)->newObj();
			record = nppr;
			break;
		}


		default:
			break;
		}
	return record;
}

void RecordMgr::deleteRecord(const Record *record)
{
	switch(_recordType) {
		case FileRecordTypeChecker::EMPTY_RECORD_TYPE:
		{
			((FreeList<EmptyRecord> *)_freeList)->deleteObj(static_cast<const EmptyRecord *>(record));
			break;
		}

		case FileRecordTypeChecker::BED3_RECORD_TYPE:
		{
			((FreeList<Bed3Interval> *)_freeList)->deleteObj(static_cast<const Bed3Interval *>(record));
			break;
		}

		case FileRecordTypeChecker::BED4_RECORD_TYPE:
		{
			((FreeList<Bed4Interval> *)_freeList)->deleteObj(static_cast<const Bed4Interval *>(record));
			break;
		}

		case FileRecordTypeChecker::BED5_RECORD_TYPE:
		{
			((FreeList<Bed5Interval> *)_freeList)->deleteObj(static_cast<const Bed5Interval *>(record));
			break;
		}

		case FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE:
		{
			((FreeList<BedGraphInterval> *)_freeList)->deleteObj(static_cast<const BedGraphInterval *>(record));
			break;
		}

		case FileRecordTypeChecker::BED6_RECORD_TYPE:
		{
			((FreeList<Bed6Interval> *)_freeList)->deleteObj(static_cast<const Bed6Interval *>(record));
			break;
		}
		case FileRecordTypeChecker::BED_PLUS_RECORD_TYPE:
		case FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE:
		{
			((FreeList<BedPlusInterval> *)_freeList)->deleteObj(static_cast<const BedPlusInterval *>(record));
			break;
		}
		case FileRecordTypeChecker::BED12_RECORD_TYPE:
		{
			((FreeList<Bed12Interval> *)_freeList)->deleteObj(static_cast<const Bed12Interval *>(record));
			break;
		}
		case FileRecordTypeChecker::BAM_RECORD_TYPE:
		{
			((FreeList<BamRecord> *)_freeList)->deleteObj(static_cast<const BamRecord *>(record));
			break;
		}

		case FileRecordTypeChecker::VCF_RECORD_TYPE:
		{
			((FreeList<VcfRecord> *)_freeList)->deleteObj(static_cast<const VcfRecord *>(record));
			break;
		}
		case FileRecordTypeChecker::GFF_RECORD_TYPE:
		{
			((FreeList<GffRecord> *)_freeList)->deleteObj(static_cast<const GffRecord *>(record));
			break;
		}
		case FileRecordTypeChecker::GFF_PLUS_RECORD_TYPE:
		{
			((FreeList<GffPlusRecord> *)_freeList)->deleteObj(static_cast<const GffPlusRecord *>(record));
			break;
		}
		case FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE:
		{
			((FreeList<NoPosPlusRecord> *)_freeList)->deleteObj(static_cast<const NoPosPlusRecord *>(record));
			break;
		}

		default:
			break;
	}
}
