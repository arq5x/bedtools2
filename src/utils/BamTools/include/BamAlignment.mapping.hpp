
struct _RefID_t {
    operator int32_t() const {return (int32_t)(_ptr()->core.tid);}
    const int32_t& operator=(const int32_t& val) {
        _ptr()->core.tid = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->RefID)))->HtsObj2();
        }
} RefID;

struct _Position_t {
    operator int32_t() const {return (int32_t)(_ptr()->core.pos);}
    const int32_t& operator=(const int32_t& val) {
        _ptr()->core.pos = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->Position)))->HtsObj2();
        }
} Position;

struct _MatePosition_t {
    operator int32_t() const {return (int32_t)(_ptr()->core.mpos);}
    const int32_t& operator=(const int32_t& val) {
        _ptr()->core.mpos = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->MatePosition)))->HtsObj2();
        }
} MatePosition;

struct _MateRefID_t {
    operator int32_t() const {return (int32_t)(_ptr()->core.mtid);}
    const int32_t& operator=(const int32_t& val) {
        _ptr()->core.mtid = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->MateRefID)))->HtsObj2();
        }
} MateRefID;

struct _InsertSize_t {
    operator int32_t() const {return (int32_t)(_ptr()->core.isize);}
    const int32_t& operator=(const int32_t& val) {
        _ptr()->core.isize = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->InsertSize)))->HtsObj2();
        }
} InsertSize;

struct _Length_t {
    operator int32_t() const {return (int32_t)(_ptr()->core.l_qseq);}
    const int32_t& operator=(const int32_t& val) {
        _ptr()->core.l_qseq = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->Length)))->HtsObj2();
        }
} Length;

struct _Bin_t {
    operator uint32_t() const {return (uint32_t)(_ptr()->core.bin);}
    const uint32_t& operator=(const uint32_t& val) {
        _ptr()->core.bin = (uint32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->Bin)))->HtsObj2();
        }
} Bin;

struct _MapQuality_t {
    operator int16_t() const {return (int16_t)(_ptr()->core.qual);}
    const int16_t& operator=(const int16_t& val) {
        _ptr()->core.qual = (int16_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->MapQuality)))->HtsObj2();
        }
} MapQuality;

struct _AlignmentFlag_t {
    operator int32_t() const {return (int32_t)(_ptr()->core.flag);}
    const int32_t& operator=(const int32_t& val) {
        _ptr()->core.flag = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->AlignmentFlag)))->HtsObj2();
        }
} AlignmentFlag;

struct _NumCigarOperations_t {
    operator uint32_t() const {return (uint32_t)(_ptr()->core.n_cigar);}
    const uint32_t& operator=(const uint32_t& val) {
        _ptr()->core.n_cigar = (uint32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->NumCigarOperations)))->HtsObj2();
        }
} NumCigarOperations;

struct _QueryNameLength_t {
    operator uint32_t() const {return (uint32_t)(_ptr()->core.l_qname);}
    const uint32_t& operator=(const uint32_t& val) {
        _ptr()->core.l_qname = (uint32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->QueryNameLength)))->HtsObj2();
        }
} QueryNameLength;

struct _QuerySequenceLength_t {
    operator uint32_t() const {return (uint32_t)(_ptr()->core.l_qseq);}
    const uint32_t& operator=(const uint32_t& val) {
        _ptr()->core.l_qseq = (uint32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr() const {
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->QuerySequenceLength)))->HtsObj2();
        }
} QuerySequenceLength;
