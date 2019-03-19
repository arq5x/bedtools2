# Makefile for htslib, a C library for high-throughput sequencing data formats.
#
#    Copyright (C) 2013-2017 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

CC     = gcc
AR     = ar
RANLIB = ranlib

# Default libraries to link if configure is not used
htslib_default_libs = -lz -lm -lbz2 -llzma

CPPFLAGS =
# TODO: probably update cram code to make it compile cleanly with -Wc++-compat
# For testing strict C99 support add -std=c99 -D_XOPEN_SOURCE=600
#CFLAGS   = -g -Wall -O2 -pedantic -std=c99 -D_XOPEN_SOURCE=600 -D__FUNCTION__=__func__
CFLAGS   = -g -Wall -O2
EXTRA_CFLAGS_PIC = -fpic
LDFLAGS  =
LIBS     = $(htslib_default_libs)

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
includedir  = $(prefix)/include
libdir      = $(exec_prefix)/lib
libexecdir  = $(exec_prefix)/libexec
datarootdir = $(prefix)/share
mandir      = $(datarootdir)/man
man1dir     = $(mandir)/man1
man5dir     = $(mandir)/man5
pkgconfigdir= $(libdir)/pkgconfig

MKDIR_P = mkdir -p
INSTALL = install -p
INSTALL_DATA    = $(INSTALL) -m 644
INSTALL_DIR     = $(MKDIR_P) -m 755
INSTALL_LIB     = $(INSTALL_DATA)
INSTALL_MAN     = $(INSTALL_DATA)
INSTALL_PROGRAM = $(INSTALL)

# Set by config.mk if plugins are enabled
plugindir =

BUILT_PROGRAMS = \
	bgzip \
	htsfile \
	tabix

BUILT_TEST_PROGRAMS = \
	test/hts_endian \
	test/fieldarith \
	test/hfile \
	test/sam \
	test/test_bgzf \
	test/test_realn \
	test/test-regidx \
	test/test_view \
	test/test-vcf-api \
	test/test-vcf-sweep \
	test/test-bcf-sr \
	test/test-bcf-translate

BUILT_THRASH_PROGRAMS = \
	test/thrash_threads1 \
	test/thrash_threads2 \
	test/thrash_threads3 \
	test/thrash_threads4 \
	test/thrash_threads5 \
	test/thrash_threads6

all: lib-static lib-shared $(BUILT_PROGRAMS) plugins $(BUILT_TEST_PROGRAMS)

HTSPREFIX =
include htslib_vars.mk

# If not using GNU make, you need to copy the version number from version.sh
# into here.
PACKAGE_VERSION := $(shell ./version.sh)
LIBHTS_SOVERSION = 2

# $(NUMERIC_VERSION) is for items that must have a numeric X.Y.Z string
# even if this is a dirty or untagged Git working tree.
NUMERIC_VERSION := $(shell ./version.sh numeric)

# Force version.h to be remade if $(PACKAGE_VERSION) has changed.
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))

version.h:
	echo '#define HTS_VERSION "$(PACKAGE_VERSION)"' > $@

print-version:
	@echo $(PACKAGE_VERSION)

show-version:
	@echo PACKAGE_VERSION = $(PACKAGE_VERSION)
	@echo NUMERIC_VERSION = $(NUMERIC_VERSION)

.SUFFIXES: .bundle .c .cygdll .dll .o .pico .so

.c.o:
	$(CC) $(CFLAGS) -I. $(CPPFLAGS) -c -o $@ $<

.c.pico:
	$(CC) $(CFLAGS) -I. $(CPPFLAGS) $(EXTRA_CFLAGS_PIC) -c -o $@ $<


LIBHTS_OBJS = \
	kfunc.o \
	knetfile.o \
	kstring.o \
	bcf_sr_sort.o \
	bgzf.o \
	errmod.o \
	faidx.o \
	hfile.o \
	hfile_net.o \
	hts.o \
	hts_os.o\
	md5.o \
	multipart.o \
	probaln.o \
	realn.o \
	regidx.o \
	sam.o \
	synced_bcf_reader.o \
	vcf_sweep.o \
	tbx.o \
	textutils.o \
	thread_pool.o \
	vcf.o \
	vcfutils.o \
	cram/cram_codecs.o \
	cram/cram_decode.o \
	cram/cram_encode.o \
	cram/cram_external.o \
	cram/cram_index.o \
	cram/cram_io.o \
	cram/cram_samtools.o \
	cram/cram_stats.o \
	cram/files.o \
	cram/mFILE.o \
	cram/open_trace_file.o \
	cram/pooled_alloc.o \
	cram/rANS_static.o \
	cram/sam_header.o \
	cram/string_alloc.o

PLUGIN_EXT  =
PLUGIN_OBJS =

cram_h = cram/cram.h $(cram_samtools_h) $(cram_sam_header_h) $(cram_structs_h) $(cram_io_h) cram/cram_encode.h cram/cram_decode.h cram/cram_stats.h cram/cram_codecs.h cram/cram_index.h $(htslib_cram_h)
cram_io_h = cram/cram_io.h $(cram_misc_h)
cram_misc_h = cram/misc.h $(cram_os_h)
cram_os_h = cram/os.h $(htslib_hts_endian_h)
cram_sam_header_h = cram/sam_header.h cram/string_alloc.h cram/pooled_alloc.h $(htslib_khash_h) $(htslib_kstring_h)
cram_samtools_h = cram/cram_samtools.h $(htslib_sam_h) $(cram_sam_header_h)
cram_structs_h = cram/cram_structs.h $(htslib_thread_pool_h) cram/string_alloc.h cram/mFILE.h $(htslib_khash_h)
cram_open_trace_file_h = cram/open_trace_file.h cram/mFILE.h
bcf_sr_sort_h = bcf_sr_sort.h $(htslib_synced_bcf_reader_h) $(htslib_kbitset_h)
hfile_internal_h = hfile_internal.h $(htslib_hfile_h) $(textutils_internal_h)
hts_internal_h = hts_internal.h $(htslib_hts_h) $(textutils_internal_h)
textutils_internal_h = textutils_internal.h $(htslib_kstring_h)
thread_pool_internal_h = thread_pool_internal.h $(htslib_thread_pool_h)


# To be effective, config.mk needs to appear after most Makefile variables are
# set but before most rules appear, so that it can both use previously-set
# variables in its own rules' prerequisites and also update variables for use
# in later rules' prerequisites.

# If your make doesn't accept -include, change this to 'include' if you are
# using the configure script or just comment the line out if you are not.
-include config.mk

# Usually config.h is generated by running configure or config.status,
# but if those aren't used create a default config.h here.
config.h:
	echo '/* Default config.h generated by Makefile */' > $@
	echo '#define HAVE_LIBBZ2 1' >> $@
	echo '#define HAVE_LIBLZMA 1' >> $@
	echo '#define HAVE_LZMA_H 1' >> $@
	echo '#define HAVE_FSEEKO 1' >> $@
	echo '#define HAVE_DRAND48 1' >> $@

# And similarly for htslib.pc.tmp ("pkg-config template").  No dependency
# on htslib.pc.in listed, as if that file is newer the usual way to regenerate
# this target is via configure or config.status rather than this rule.
htslib.pc.tmp:
	sed -e '/^static_libs=/s/@static_LIBS@/$(htslib_default_libs)/;s#@[^-][^@]*@##g' htslib.pc.in > $@

# Create a makefile fragment listing the libraries and LDFLAGS needed for
# static linking.  This can be included by projects that want to build
# and link against the htslib source tree instead of an installed library.
htslib_static.mk: htslib.pc.tmp
	sed -n '/^static_libs=/s/[^=]*=/HTSLIB_static_LIBS = /p;/^static_ldflags=/s/[^=]*=/HTSLIB_static_LDFLAGS = /p' $< > $@


lib-static: libhts.a

# $(shell), :=, and ifeq/.../endif are GNU Make-specific.  If you don't have
# GNU Make, comment out the parts of these conditionals that don't apply.
ifneq "$(origin PLATFORM)" "file"
PLATFORM := $(shell uname -s)
endif
ifeq "$(PLATFORM)" "Darwin"
SHLIB_FLAVOUR = dylib
lib-shared: libhts.dylib
else ifeq "$(findstring CYGWIN,$(PLATFORM))" "CYGWIN"
SHLIB_FLAVOUR = cygdll
lib-shared: cyghts-$(LIBHTS_SOVERSION).dll
else ifeq "$(findstring MSYS,$(PLATFORM))" "MSYS"
SHLIB_FLAVOUR = dll
lib-shared: hts-$(LIBHTS_SOVERSION).dll
else
SHLIB_FLAVOUR = so
lib-shared: libhts.so
endif

BUILT_PLUGINS = $(PLUGIN_OBJS:.o=$(PLUGIN_EXT))

plugins: $(BUILT_PLUGINS)


libhts.a: $(LIBHTS_OBJS)
	@-rm -f $@
	$(AR) -rc $@ $(LIBHTS_OBJS)
	-$(RANLIB) $@

print-config:
	@echo LDFLAGS = $(LDFLAGS)
	@echo LIBHTS_OBJS = $(LIBHTS_OBJS)
	@echo LIBS = $(LIBS)
	@echo PLATFORM = $(PLATFORM)

# The target here is libhts.so, as that is the built file that other rules
# depend upon and that is used when -lhts appears in other program's recipes.
# As a byproduct invisible to make, libhts.so.NN is also created, as it is the
# file used at runtime (when $LD_LIBRARY_PATH includes the build directory).

libhts.so: $(LIBHTS_OBJS:.o=.pico)
	$(CC) -shared -Wl,-soname,libhts.so.$(LIBHTS_SOVERSION) $(LDFLAGS) -o $@ $(LIBHTS_OBJS:.o=.pico) $(LIBS) -lpthread
	ln -sf $@ libhts.so.$(LIBHTS_SOVERSION)

# Similarly this also creates libhts.NN.dylib as a byproduct, so that programs
# when run can find this uninstalled shared library (when $DYLD_LIBRARY_PATH
# includes this project's build directory).

libhts.dylib: $(LIBHTS_OBJS)
	$(CC) -dynamiclib -install_name $(libdir)/libhts.$(LIBHTS_SOVERSION).dylib -current_version $(NUMERIC_VERSION) -compatibility_version $(LIBHTS_SOVERSION) $(LDFLAGS) -o $@ $(LIBHTS_OBJS) $(LIBS)
	ln -sf $@ libhts.$(LIBHTS_SOVERSION).dylib

cyghts-$(LIBHTS_SOVERSION).dll: $(LIBHTS_OBJS)
	$(CC) -shared -Wl,--out-implib=libhts.dll.a -Wl,--export-all-symbols -Wl,--enable-auto-import $(LDFLAGS) -o $@ -Wl,--whole-archive $(LIBHTS_OBJS) -Wl,--no-whole-archive $(LIBS) -lpthread

hts-$(LIBHTS_SOVERSION).dll: $(LIBHTS_OBJS)
	$(CC) -shared -Wl,--out-implib=hts.dll.a -Wl,--export-all-symbols -Wl,--enable-auto-import $(LDFLAGS) -o $@ -Wl,--whole-archive $(LIBHTS_OBJS) -Wl,--no-whole-archive $(LIBS) -lpthread


.pico.so:
	$(CC) -shared -Wl,-E $(LDFLAGS) -o $@ $< $(LIBS) -lpthread

.o.bundle:
	$(CC) -bundle -Wl,-undefined,dynamic_lookup $(LDFLAGS) -o $@ $< $(LIBS)

.o.cygdll:
	$(CC) -shared $(LDFLAGS) -o $@ $< libhts.dll.a $(LIBS)

.o.dll:
	$(CC) -shared $(LDFLAGS) -o $@ $< hts.dll.a $(LIBS)


bgzf.o bgzf.pico: bgzf.c config.h $(htslib_hts_h) $(htslib_bgzf_h) $(htslib_hfile_h) $(htslib_thread_pool_h) $(htslib_hts_endian_h) cram/pooled_alloc.h $(htslib_khash_h)
errmod.o errmod.pico: errmod.c config.h $(htslib_hts_h) $(htslib_ksort_h) $(htslib_hts_os_h)
kstring.o kstring.pico: kstring.c config.h $(htslib_kstring_h)
knetfile.o knetfile.pico: knetfile.c config.h $(htslib_hts_log_h) $(htslib_knetfile_h)
hfile.o hfile.pico: hfile.c config.h $(htslib_hfile_h) $(hfile_internal_h) $(hts_internal_h) $(htslib_khash_h)
hfile_gcs.o hfile_gcs.pico: hfile_gcs.c config.h $(htslib_hts_h) $(htslib_kstring_h) $(hfile_internal_h)
hfile_libcurl.o hfile_libcurl.pico: hfile_libcurl.c config.h $(hfile_internal_h) $(htslib_hts_h) $(htslib_kstring_h) $(htslib_khash_h)
hfile_net.o hfile_net.pico: hfile_net.c config.h $(hfile_internal_h) $(htslib_knetfile_h)
hfile_s3.o hfile_s3.pico: hfile_s3.c config.h $(hfile_internal_h) $(htslib_hts_h) $(htslib_kstring_h)
hts.o hts.pico: hts.c config.h $(htslib_hts_h) $(htslib_bgzf_h) $(cram_h) $(htslib_hfile_h) $(htslib_hts_endian_h) version.h $(hts_internal_h) $(hfile_internal_h) $(htslib_hts_os_h) $(htslib_khash_h) $(htslib_kseq_h) $(htslib_ksort_h)
hts_os.o hts_os.pico: hts_os.c config.h os/rand.c
vcf.o vcf.pico: vcf.c config.h $(htslib_vcf_h) $(htslib_bgzf_h) $(htslib_tbx_h) $(htslib_hfile_h) $(hts_internal_h) $(htslib_khash_str2int_h) $(htslib_kstring_h) $(htslib_khash_h) $(htslib_kseq_h) $(htslib_hts_endian_h)
sam.o sam.pico: sam.c config.h $(htslib_sam_h) $(htslib_bgzf_h) $(cram_h) $(hts_internal_h) $(htslib_hfile_h) $(htslib_khash_h) $(htslib_kseq_h) $(htslib_kstring_h) $(htslib_hts_endian_h)
tbx.o tbx.pico: tbx.c config.h $(htslib_tbx_h) $(htslib_bgzf_h) $(htslib_hts_endian_h) $(hts_internal_h) $(htslib_khash_h)
faidx.o faidx.pico: faidx.c config.h $(htslib_bgzf_h) $(htslib_faidx_h) $(htslib_hfile_h) $(htslib_khash_h) $(htslib_kstring_h) $(hts_internal_h)
bcf_sr_sort.o bcf_sr_sort.pico: bcf_sr_sort.c config.h $(bcf_sr_sort_h) $(htslib_khash_str2int_h)
synced_bcf_reader.o synced_bcf_reader.pico: synced_bcf_reader.c config.h $(htslib_synced_bcf_reader_h) $(htslib_kseq_h) $(htslib_khash_str2int_h) $(htslib_bgzf_h) $(htslib_thread_pool_h) $(bcf_sr_sort_h)
vcf_sweep.o vcf_sweep.pico: vcf_sweep.c config.h $(htslib_vcf_sweep_h) $(htslib_bgzf_h)
vcfutils.o vcfutils.pico: vcfutils.c config.h $(htslib_vcfutils_h) $(htslib_kbitset_h)
kfunc.o kfunc.pico: kfunc.c config.h $(htslib_kfunc_h)
regidx.o regidx.pico: regidx.c config.h $(htslib_hts_h) $(htslib_kstring_h) $(htslib_kseq_h) $(htslib_khash_str2int_h) $(htslib_regidx_h) $(hts_internal_h)
md5.o md5.pico: md5.c config.h $(htslib_hts_h) $(htslib_hts_endian_h)
multipart.o multipart.pico: multipart.c config.h $(htslib_kstring_h) $(hts_internal_h) $(hfile_internal_h)
plugin.o plugin.pico: plugin.c config.h $(hts_internal_h) $(htslib_kstring_h)
probaln.o probaln.pico: probaln.c config.h $(htslib_hts_h)
realn.o realn.pico: realn.c config.h $(htslib_hts_h) $(htslib_sam_h)
textutils.o textutils.pico: textutils.c config.h $(htslib_hfile_h) $(htslib_kstring_h) $(hts_internal_h)

cram/cram_codecs.o cram/cram_codecs.pico: cram/cram_codecs.c config.h $(cram_h)
cram/cram_decode.o cram/cram_decode.pico: cram/cram_decode.c config.h $(cram_h) $(cram_os_h) $(htslib_hts_h)
cram/cram_encode.o cram/cram_encode.pico: cram/cram_encode.c config.h $(cram_h) $(cram_os_h) $(htslib_hts_h) $(htslib_hts_endian_h)
cram/cram_external.o cram/cram_external.pico: cram/cram_external.c config.h $(htslib_hfile_h) $(cram_h)
cram/cram_index.o cram/cram_index.pico: cram/cram_index.c config.h $(htslib_bgzf_h) $(htslib_hfile_h) $(hts_internal_h) $(cram_h) $(cram_os_h)
cram/cram_io.o cram/cram_io.pico: cram/cram_io.c config.h os/lzma_stub.h $(cram_h) $(cram_os_h) $(htslib_hts_h) $(cram_open_trace_file_h) cram/rANS_static.h $(htslib_hfile_h) $(htslib_bgzf_h) $(htslib_faidx_h) $(hts_internal_h)
cram/cram_samtools.o cram/cram_samtools.pico: cram/cram_samtools.c config.h $(cram_h) $(htslib_sam_h)
cram/cram_stats.o cram/cram_stats.pico: cram/cram_stats.c config.h $(cram_h) $(cram_os_h)
cram/files.o cram/files.pico: cram/files.c config.h $(cram_misc_h)
cram/mFILE.o cram/mFILE.pico: cram/mFILE.c config.h $(htslib_hts_log_h) $(cram_os_h) cram/mFILE.h
cram/open_trace_file.o cram/open_trace_file.pico: cram/open_trace_file.c config.h $(cram_os_h) $(cram_open_trace_file_h) $(cram_misc_h) $(htslib_hfile_h) $(htslib_hts_log_h)
cram/pooled_alloc.o cram/pooled_alloc.pico: cram/pooled_alloc.c config.h cram/pooled_alloc.h $(cram_misc_h)
cram/rANS_static.o cram/rANS_static.pico: cram/rANS_static.c config.h cram/rANS_static.h cram/rANS_byte.h
cram/sam_header.o cram/sam_header.pico: cram/sam_header.c config.h $(htslib_hts_log_h) $(cram_sam_header_h) cram/string_alloc.h
cram/string_alloc.o cram/string_alloc.pico: cram/string_alloc.c config.h cram/string_alloc.h
thread_pool.o thread_pool.pico: thread_pool.c config.h $(thread_pool_internal_h)


bgzip: bgzip.o libhts.a
	$(CC) $(LDFLAGS) -o $@ bgzip.o libhts.a $(LIBS) -lpthread

htsfile: htsfile.o libhts.a
	$(CC) $(LDFLAGS) -o $@ htsfile.o libhts.a $(LIBS) -lpthread

tabix: tabix.o libhts.a
	$(CC) $(LDFLAGS) -o $@ tabix.o libhts.a $(LIBS) -lpthread

bgzip.o: bgzip.c config.h $(htslib_bgzf_h) $(htslib_hts_h)
htsfile.o: htsfile.c config.h $(htslib_hfile_h) $(htslib_hts_h) $(htslib_sam_h) $(htslib_vcf_h)
tabix.o: tabix.c config.h $(htslib_tbx_h) $(htslib_sam_h) $(htslib_vcf_h) $(htslib_kseq_h) $(htslib_bgzf_h) $(htslib_hts_h) $(htslib_regidx_h)


# For tests that might use it, set $REF_PATH explicitly to use only reference
# areas within the test suite (or set it to ':' to use no reference areas).
#
# If using MSYS, avoid poor shell expansion via:
#    MSYS2_ARG_CONV_EXCL="*" make check
check test: $(BUILT_PROGRAMS) $(BUILT_TEST_PROGRAMS)
	test/hts_endian
	test/fieldarith test/fieldarith.sam
	test/hfile
	test/test_bgzf test/bgziptest.txt
	cd test/tabix && ./test-tabix.sh tabix.tst
	REF_PATH=: test/sam test/ce.fa test/faidx.fa test/fastqs.fq
	test/test-regidx
	cd test && REF_PATH=: ./test.pl $${TEST_OPTS:-}

test/hts_endian: test/hts_endian.o
	$(CC) $(LDFLAGS) -o $@ test/hts_endian.o $(LIBS)

test/fieldarith: test/fieldarith.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/fieldarith.o libhts.a $(LIBS) -lpthread

test/hfile: test/hfile.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/hfile.o libhts.a $(LIBS) -lpthread

test/sam: test/sam.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/sam.o libhts.a $(LIBS) -lpthread

test/test_bgzf: test/test_bgzf.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/test_bgzf.o libhts.a -lz $(LIBS) -lpthread

test/test_realn: test/test_realn.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/test_realn.o libhts.a $(LIBS) -lpthread

test/test-regidx: test/test-regidx.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/test-regidx.o libhts.a $(LIBS) -lpthread

test/test_view: test/test_view.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/test_view.o libhts.a $(LIBS) -lpthread

test/test-vcf-api: test/test-vcf-api.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/test-vcf-api.o libhts.a $(LIBS) -lpthread

test/test-vcf-sweep: test/test-vcf-sweep.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/test-vcf-sweep.o libhts.a $(LIBS) -lpthread

test/test-bcf-sr: test/test-bcf-sr.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/test-bcf-sr.o libhts.a -lz $(LIBS) -lpthread

test/test-bcf-translate: test/test-bcf-translate.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/test-bcf-translate.o libhts.a -lz $(LIBS) -lpthread

test/hts_endian.o: test/hts_endian.c config.h $(htslib_hts_endian_h)
test/fieldarith.o: test/fieldarith.c config.h $(htslib_sam_h)
test/hfile.o: test/hfile.c config.h $(htslib_hfile_h) $(htslib_hts_defs_h)
test/sam.o: test/sam.c config.h $(htslib_hts_defs_h) $(htslib_sam_h) $(htslib_faidx_h) $(htslib_kstring_h)
test/test_bgzf.o: test/test_bgzf.c config.h $(htslib_bgzf_h) $(htslib_hfile_h) $(hfile_internal_h)
test/test-realn.o: test/test_realn.c config.h $(htslib_hts_h) $(htslib_sam_h) $(htslib_faidx_h)
test/test-regidx.o: test/test-regidx.c config.h $(htslib_regidx_h) $(hts_internal_h)
test/test_view.o: test/test_view.c config.h $(cram_h) $(htslib_sam_h)
test/test-vcf-api.o: test/test-vcf-api.c config.h $(htslib_hts_h) $(htslib_vcf_h) $(htslib_kstring_h) $(htslib_kseq_h)
test/test-vcf-sweep.o: test/test-vcf-sweep.c config.h $(htslib_vcf_sweep_h)
test/test-bcf-sr.o: test/test-bcf-sr.c config.h $(htslib_synced_bcf_reader_h)
test/test-bcf-translate.o: test/test-bcf-translate.c config.h $(htslib_vcf_h)


test/thrash_threads1: test/thrash_threads1.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/thrash_threads1.o libhts.a -lz $(LIBS) -lpthread

test/thrash_threads2: test/thrash_threads2.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/thrash_threads2.o libhts.a -lz $(LIBS) -lpthread

test/thrash_threads3: test/thrash_threads3.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/thrash_threads3.o libhts.a -lz $(LIBS) -lpthread

test/thrash_threads4: test/thrash_threads4.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/thrash_threads4.o libhts.a -lz $(LIBS) -lpthread

test/thrash_threads5: test/thrash_threads5.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/thrash_threads5.o libhts.a -lz $(LIBS) -lpthread

test/thrash_threads6: test/thrash_threads6.o libhts.a
	$(CC) $(LDFLAGS) -o $@ test/thrash_threads6.o libhts.a -lz $(LIBS) -lpthread

test_thrash: $(BUILT_THRASH_PROGRAMS)


install: libhts.a $(BUILT_PROGRAMS) $(BUILT_PLUGINS) installdirs install-$(SHLIB_FLAVOUR) install-pkgconfig
	$(INSTALL_PROGRAM) $(BUILT_PROGRAMS) $(DESTDIR)$(bindir)
	if test -n "$(BUILT_PLUGINS)"; then $(INSTALL_PROGRAM) $(BUILT_PLUGINS) $(DESTDIR)$(plugindir); fi
	$(INSTALL_DATA) htslib/*.h $(DESTDIR)$(includedir)/htslib
	$(INSTALL_DATA) libhts.a $(DESTDIR)$(libdir)/libhts.a
	$(INSTALL_MAN) bgzip.1 htsfile.1 tabix.1 $(DESTDIR)$(man1dir)
	$(INSTALL_MAN) faidx.5 sam.5 vcf.5 $(DESTDIR)$(man5dir)

installdirs:
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(includedir) $(DESTDIR)$(includedir)/htslib $(DESTDIR)$(libdir) $(DESTDIR)$(man1dir) $(DESTDIR)$(man5dir) $(DESTDIR)$(pkgconfigdir)
	if test -n "$(plugindir)"; then $(INSTALL_DIR) $(DESTDIR)$(plugindir); fi

# After installation, the real file in $(libdir) will be libhts.so.X.Y.Z,
# with symlinks libhts.so (used via -lhts during linking of client programs)
# and libhts.so.NN (used by client executables at runtime).

install-so: libhts.so installdirs
	$(INSTALL_LIB) libhts.so $(DESTDIR)$(libdir)/libhts.so.$(PACKAGE_VERSION)
	ln -sf libhts.so.$(PACKAGE_VERSION) $(DESTDIR)$(libdir)/libhts.so
	ln -sf libhts.so.$(PACKAGE_VERSION) $(DESTDIR)$(libdir)/libhts.so.$(LIBHTS_SOVERSION)

install-cygdll: cyghts-$(LIBHTS_SOVERSION).dll installdirs
	$(INSTALL_PROGRAM) cyghts-$(LIBHTS_SOVERSION).dll $(DESTDIR)$(bindir)/cyghts-$(LIBHTS_SOVERSION).dll
	$(INSTALL_PROGRAM) libhts.dll.a $(DESTDIR)$(libdir)/libhts.dll.a

install-dll: hts-$(LIBHTS_SOVERSION).dll installdirs
	$(INSTALL_PROGRAM) hts-$(LIBHTS_SOVERSION).dll $(DESTDIR)$(bindir)/hts-$(LIBHTS_SOVERSION).dll
	$(INSTALL_PROGRAM) hts.dll.a $(DESTDIR)$(libdir)/hts.dll.a

install-dylib: libhts.dylib installdirs
	$(INSTALL_PROGRAM) libhts.dylib $(DESTDIR)$(libdir)/libhts.$(PACKAGE_VERSION).dylib
	ln -sf libhts.$(PACKAGE_VERSION).dylib $(DESTDIR)$(libdir)/libhts.dylib
	ln -sf libhts.$(PACKAGE_VERSION).dylib $(DESTDIR)$(libdir)/libhts.$(LIBHTS_SOVERSION).dylib

# Substitute these pseudo-autoconf variables only at install time
# so that "make install prefix=/prefix/path" etc continue to work.
install-pkgconfig: htslib.pc.tmp installdirs
	sed -e 's#@-includedir@#$(includedir)#g;s#@-libdir@#$(libdir)#g;s#@-PACKAGE_VERSION@#$(PACKAGE_VERSION)#g' htslib.pc.tmp > $(DESTDIR)$(pkgconfigdir)/htslib.pc
	chmod 644 $(DESTDIR)$(pkgconfigdir)/htslib.pc

# A pkg-config file (suitable for copying to $PKG_CONFIG_PATH) that provides
# flags for building against the uninstalled library in this build directory.
htslib-uninstalled.pc: htslib.pc.tmp
	sed -e 's#@-includedir@#'`pwd`'#g;s#@-libdir@#'`pwd`'#g' htslib.pc.tmp > $@


testclean:
	-rm -f test/*.tmp test/*.tmp.* test/tabix/*.tmp.* test/tabix/FAIL*

mostlyclean: testclean
	-rm -f *.o *.pico cram/*.o cram/*.pico test/*.o test/*.dSYM version.h

clean: mostlyclean clean-$(SHLIB_FLAVOUR)
	-rm -f libhts.a $(BUILT_PROGRAMS) $(BUILT_PLUGINS) $(BUILT_TEST_PROGRAMS) $(BUILT_THRASH_PROGRAMS)

distclean maintainer-clean: clean
	-rm -f config.cache config.h config.log config.mk config.status
	-rm -f TAGS *.pc.tmp *-uninstalled.pc htslib_static.mk
	-rm -rf autom4te.cache

clean-so:
	-rm -f libhts.so libhts.so.*

clean-cygdll:
	-rm -f cyghts-*.dll libhts.dll.a

clean-dll:
	-rm -f hts-*.dll hts.dll.a

clean-dylib:
	-rm -f libhts.dylib libhts.*.dylib


tags TAGS:
	ctags -f TAGS *.[ch] cram/*.[ch] htslib/*.h

# We recommend libhts-using programs be built against a separate htslib
# installation.  However if you feel that you must bundle htslib source
# code with your program, this hook enables Automake-style "make dist"
# for this subdirectory.  If you do bundle an htslib snapshot, please
# add identifying information to $(PACKAGE_VERSION) as appropriate.
# (The wildcards attempt to omit non-exported files (.git*, README.md,
# etc) and other detritus that might be in the top-level directory.)
distdir:
	@if [ -z "$(distdir)" ]; then echo "Please supply a distdir=DIR argument."; false; fi
	tar -c *.[ch15] [ILMNRchtv]*[ELSbcekmnth] | (cd $(distdir) && tar -x)
	+cd $(distdir) && $(MAKE) distclean

force:


.PHONY: all check clean distclean distdir force
.PHONY: install install-pkgconfig installdirs lib-shared lib-static
.PHONY: maintainer-clean mostlyclean plugins print-config print-version
.PHONY: show-version tags test testclean
.PHONY: clean-so install-so
.PHONY: clean-cygdll install-cygdll
.PHONY: clean-dll install-dll
.PHONY: clean-dylib install-dylib
