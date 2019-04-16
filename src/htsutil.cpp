#include <iostream>
#include <string>

#include <htslib/sam.h>

int convert(htsExactFormat outformat, const char *outfname, bool outheaders,
            htsExactFormat informat, const char *infname,
            const char *reffname = NULL)
{
    htsFile* in = NULL;
    htsFile* out = NULL;
    bam_hdr_t* hdr = NULL;
    bam1_t* b = NULL;

    try {
        in = hts_open(infname, "r");
        if (in == NULL)
            throw std::string("can't open ") + infname;

        if (hts_get_format(in)->format != informat)
            throw std::string(infname) + ": file is in the wrong format";

        const char *outf = "w";
        switch (outformat) {
        case sam:  outf = "w"; break;
        case bam:  outf = "wb"; break;
        case cram: outf = "wc"; break;
        default: throw std::string("invalid outformat argument");
        }

        out = hts_open(outfname, outf);
        if (out == NULL)
            throw std::string("can't write to ") + outfname;

        if (reffname) {
            hts_set_fai_filename(in, reffname);
            hts_set_fai_filename(out, reffname);
        }

        hdr = sam_hdr_read(in);
        if (hdr == NULL)
            throw std::string("can't read headers from ") + infname;

        if (outheaders) {
            if (sam_hdr_write(out, hdr) < 0)
                throw std::string("can't write headers to ") + outfname;
        }

        int ret;
        b = bam_init1();
        while ((ret = sam_read1(in, hdr, b)) >= 0) {
            if (sam_write1(out, hdr, b) < 0)
                throw std::string("can't write record to ") + outfname;
        }
        if (ret < -1)
            throw std::string("truncated input file");

        bam_destroy1(b);
        bam_hdr_destroy(hdr);
        hts_close(out);
        hts_close(in);
    }
    catch (const std::string& what) {
        std::cerr << "htsutil: " << what << '\n';
        bam_destroy1(b);
        bam_hdr_destroy(hdr);
        if (out) hts_close(out);
        if (in) hts_close(in);
        return 1;
    }

    return 0;
}

int index(const char *fname, bool csi)
{
    if (std::string(fname) == "-") {
        std::cerr << "htsutil: can't index standard output\n";
        return 1;
    }

    int ret = sam_index_build(fname, csi? 14 : 0);
    if (ret < 0) return -ret; // Convert to exit status
    return 0;
}

int main(int argc, char **argv)
{
    std::string cmd = (argc >= 2)? argv[1] : "--help";
    const char *src = (argc >= 3)? argv[2] : "-";
    const char *dest= (argc >= 4)? argv[3] : "-";

    if (cmd == "viewbamrecords")
        return convert(sam, "-", false, bam, src);
    else if (cmd == "viewbam")
        return convert(sam, "-", true, bam, src);
    else if (cmd == "viewcramrecords")
        return convert(sam, "-", false, cram, src, (argc >= 4)? argv[3] : NULL);
    else if (cmd == "index")
        return index(src, false);
    else if (cmd == "indexcsi")
        return index(src, true);
    else if (cmd == "samtobam")
        return convert(bam, dest, true, sam, src);
    else if (cmd == "samtoindexedbam") {
        int ret = convert(bam, dest, true, sam, src);
        if (ret != 0) return ret;
        return index(dest, false);
    }
    else {
        std::cerr <<
"Usage: htsutil index FILE\n"
"       htsutil indexcsi FILE\n"
"       htsutil samtobam [SRCFILE [DESTFILE]]\n"
"       htsutil samtoindexedbam SRCFILE DESTFILE\n"
"       htsutil viewbam [FILE]\n"
"       htsutil viewbamrecords [FILE]\n"
"       htsutil viewcramrecords [FILE [REFFILE]]\n";
        return 1;
    }
}
