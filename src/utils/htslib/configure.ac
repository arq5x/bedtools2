# Configure script for htslib, a C library for high-throughput sequencing data.
#
#    Copyright (C) 2015-2017 Genome Research Ltd.
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

dnl Process this file with autoconf to produce a configure script
AC_INIT([HTSlib], m4_esyscmd_s([./version.sh 2>/dev/null]),
        [samtools-help@lists.sourceforge.net], [], [http://www.htslib.org/])
AC_PREREQ(2.63)  dnl This version introduced 4-argument AC_CHECK_HEADER
AC_CONFIG_SRCDIR(hts.c)
AC_CONFIG_HEADERS(config.h)

dnl Copyright notice to be copied into the generated configure script
AC_COPYRIGHT([Portions copyright (C) 2016 Genome Research Ltd.

This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

dnl Notes to be copied (by autoheader) into the generated config.h.in
AH_TOP([/* If you use configure, this file provides @%:@defines reflecting your
   configuration choices.  If you have not run configure, suitable
   conservative defaults will be used.

   Autoheader adds a number of items to this template file that are not
   used by HTSlib: STDC_HEADERS and most HAVE_*_H header file defines
   are immaterial, as we assume standard ISO C headers and facilities;
   the PACKAGE_* defines are unused and are overridden by the more
   accurate PACKAGE_VERSION as computed by the Makefile.  */])

AC_PROG_CC
AC_PROG_RANLIB

dnl Avoid chicken-and-egg problem where pkg-config supplies the
dnl PKG_PROG_PKG_CONFIG macro, but we want to use it to check
dnl for pkg-config...
m4_ifdef([PKG_PROG_PKG_CONFIG], [PKG_PROG_PKG_CONFIG], [PKG_CONFIG=""])

need_crypto=no
pc_requires=
static_LDFLAGS=
static_LIBS='-lz -lm'
private_LIBS=

AC_ARG_ENABLE([bz2],
  [AS_HELP_STRING([--disable-bz2],
                  [omit support for BZ2-compressed CRAM files])],
  [], [enable_bz2=yes])

AC_ARG_ENABLE([gcs],
  [AS_HELP_STRING([--enable-gcs],
                  [support Google Cloud Storage URLs])],
  [], [enable_gcs=check])

AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

AC_ARG_ENABLE([libcurl],
  [AS_HELP_STRING([--enable-libcurl],
                  [enable libcurl-based support for http/https/etc URLs])],
  [], [enable_libcurl=check])

AC_ARG_ENABLE([lzma],
  [AS_HELP_STRING([--disable-lzma],
                  [omit support for LZMA-compressed CRAM files])],
  [], [enable_lzma=yes])

AC_ARG_ENABLE([plugins],
  [AS_HELP_STRING([--enable-plugins],
                  [enable separately-compiled plugins for file access])],
  [], [enable_plugins=no])
AC_SUBST(enable_plugins)

AC_ARG_WITH([plugin-dir],
  [AS_HELP_STRING([--with-plugin-dir=DIR],
                  [plugin installation location [LIBEXECDIR/htslib]])],
  [case $withval in
     yes|no) AC_MSG_ERROR([no directory specified for --with-plugin-dir]) ;;
   esac],
   [with_plugin_dir='$(libexecdir)/htslib'])
AC_SUBST([plugindir], $with_plugin_dir)

AC_ARG_WITH([plugin-path],
  [AS_HELP_STRING([--with-plugin-path=PATH],
                  [default HTS_PATH plugin search path [PLUGINDIR]])],
  [case $withval in
     yes) AC_MSG_ERROR([no path specified for --with-plugin-path]) ;;
     no)  with_plugin_path= ;;
   esac],
  [with_plugin_path=$with_plugin_dir])
AC_SUBST([pluginpath], $with_plugin_path)

AC_ARG_ENABLE([s3],
  [AS_HELP_STRING([--enable-s3],
                  [support Amazon AWS S3 URLs])],
  [], [enable_s3=check])

test -n "$host_alias" || host_alias=unknown-`uname -s`
AC_MSG_CHECKING([shared library type for $host_alias])
case $host_alias in
  *-cygwin* | *-CYGWIN*)
    host_result="Cygwin DLL"
    PLATFORM=CYGWIN
    PLUGIN_EXT=.cygdll
    ;;
  *-darwin* | *-Darwin*)
    host_result="Darwin dylib"
    PLATFORM=Darwin
    PLUGIN_EXT=.bundle
    ;;
  *-msys* | *-MSYS* | *-mingw* | *-MINGW*)
    host_result="MSYS dll"
    PLATFORM=MSYS
    PLUGIN_EXT=.dll
    # This also sets __USE_MINGW_ANSI_STDIO which in turn makes PRId64,
    # %lld and %z printf formats work.  It also enforces the snprintf to
    # be C99 compliant so it returns the correct values (in kstring.c).
    CPPFLAGS="$CPPCFLAGS -D_XOPEN_SOURCE=600"
    ;;
  *)
    host_result="plain .so"
    PLATFORM=default
    PLUGIN_EXT=.so
    ;;
esac
AC_MSG_RESULT([$host_result])
AC_SUBST([PLATFORM])

dnl FIXME This pulls in dozens of standard header checks
AC_FUNC_MMAP
AC_CHECK_FUNCS([gmtime_r fsync drand48])

# Darwin has a dubious fdatasync() symbol, but no declaration in <unistd.h>
AC_CHECK_DECL([fdatasync(int)], [AC_CHECK_FUNCS(fdatasync)])

if test $enable_plugins != no; then
  AC_SEARCH_LIBS([dlopen], [dl], [],
    [AC_MSG_ERROR([dlopen() not found

Plugin support requires dynamic linking facilities from the operating system.
Either configure with --disable-plugins or resolve this error to build HTSlib.])])
  # TODO Test whether this is required and/or needs tweaking per-platform
  LDFLAGS="$LDFLAGS -rdynamic"
  static_LDFLAGS="$static_LDFLAGS -rdynamic"
  case "$ac_cv_search_dlopen" in
    -l*) static_LIBS="$static_LIBS $ac_cv_search_dlopen" ;;
  esac
  AC_DEFINE([ENABLE_PLUGINS], 1, [Define if HTSlib should enable plugins.])
  AC_SUBST([PLUGIN_EXT])
  AC_DEFINE_UNQUOTED([PLUGIN_EXT], ["$PLUGIN_EXT"],
                     [Platform-dependent plugin filename extension.])
fi

AC_SEARCH_LIBS([log], [m], [],
  [AC_MSG_ERROR([log() not found

HTSLIB requires a working floating-point math library.
FAILED.  This error must be resolved in order to build HTSlib successfully.])])

zlib_devel=ok
dnl Set a trivial non-empty INCLUDES to avoid excess default includes tests
AC_CHECK_HEADER([zlib.h], [], [zlib_devel=missing], [;])
AC_CHECK_LIB(z, inflate,  [], [zlib_devel=missing])

if test $zlib_devel != ok; then
  AC_MSG_ERROR([zlib development files not found

HTSlib uses compression routines from the zlib library <http://zlib.net>.
Building HTSlib requires zlib development files to be installed on the build
machine; you may need to ensure a package such as zlib1g-dev (on Debian or
Ubuntu Linux) or zlib-devel (on RPM-based Linux distributions or Cygwin)
is installed.

FAILED.  This error must be resolved in order to build HTSlib successfully.])
fi

dnl connect() etc. fns are in libc on linux, but libsocket on illumos/Solaris
AC_SEARCH_LIBS([recv], [socket ws2_32], [
if test "$ac_cv_search_recv" != "none required"
then
  static_LIBS="$static_LIBS $ac_cv_search_recv"
fi], [AC_MSG_ERROR([unable to find the recv() function])])

if test "$enable_bz2" != no; then
  bz2_devel=ok
  AC_CHECK_HEADER([bzlib.h], [], [bz2_devel=missing], [;])
  AC_CHECK_LIB([bz2], [BZ2_bzBuffToBuffCompress], [], [bz2_devel=missing])
  if test $bz2_devel != ok; then
    AC_MSG_ERROR([libbzip2 development files not found

The CRAM format may use bzip2 compression, which is implemented in HTSlib
by using compression routines from libbzip2 <http://www.bzip.org/>.

Building HTSlib requires libbzip2 development files to be installed on the
build machine; you may need to ensure a package such as libbz2-dev (on Debian
or Ubuntu Linux) or bzip2-devel (on RPM-based Linux distributions or Cygwin)
is installed.

Either configure with --disable-bz2 (which will make some CRAM files
produced elsewhere unreadable) or resolve this error to build HTSlib.])
  fi
dnl Unfortunately the 'bzip2' package-cfg module is not standard.
dnl Redhat/Fedora has it; Debian/Ubuntu does not.
  if test -n "$PKG_CONFIG" && "$PKG_CONFIG" --exists bzip2; then
     pc_requires="$pc_requires bzip2"
  else
     private_LIBS="$private_LIBS -lbz2"
  fi  
  static_LIBS="$static_LIBS -lbz2"
fi

if test "$enable_lzma" != no; then
  lzma_devel=ok
  AC_CHECK_HEADER([lzma.h], [], [lzma_devel=missing], [;])
  AC_CHECK_LIB([lzma], [lzma_easy_buffer_encode], [], [lzma_devel=missing])
  if test $lzma_devel != ok; then
    AC_MSG_ERROR([liblzma development files not found

The CRAM format may use LZMA2 compression, which is implemented in HTSlib
by using compression routines from liblzma <http://tukaani.org/xz/>.

Building HTSlib requires liblzma development files to be installed on the
build machine; you may need to ensure a package such as liblzma-dev (on Debian
or Ubuntu Linux), xz-devel (on RPM-based Linux distributions or Cygwin), or
xz (via Homebrew on macOS) is installed; or build XZ Utils from source.

Either configure with --disable-lzma (which will make some CRAM files
produced elsewhere unreadable) or resolve this error to build HTSlib.])
  fi
  pc_requires="$pc_requires liblzma"
  static_LIBS="$static_LIBS -llzma"
fi

libcurl=disabled
if test "$enable_libcurl" != no; then
  AC_CHECK_LIB([curl], [curl_easy_pause],
    [AC_DEFINE([HAVE_LIBCURL], 1, [Define if libcurl file access is enabled.])
     libcurl=enabled],
    [AC_CHECK_LIB([curl], [curl_easy_init],
       [message="library is too old (7.18+ required)"],
       [message="library not found"])
     case "$enable_libcurl" in
       check) AC_MSG_WARN([libcurl not enabled: $message]) ;;
       *) AC_MSG_ERROR([libcurl $message

Support for HTTPS and other SSL-based URLs requires routines from the libcurl
library <http://curl.haxx.se/libcurl/>.  Building HTSlib with libcurl enabled
requires libcurl development files to be installed on the build machine; you
may need to ensure a package such as libcurl4-{gnutls,nss,openssl}-dev (on
Debian or Ubuntu Linux) or libcurl-devel (on RPM-based Linux distributions
or Cygwin) is installed.

Either configure with --disable-libcurl or resolve this error to build HTSlib.])
       ;;
     esac])
dnl -lcurl is only needed for static linking if hfile_libcurl is not a plugin
  if test "$libcurl" = enabled ; then
    if test "$enable_plugins" != yes ; then
      static_LIBS="$static_LIBS -lcurl"
    fi
  fi
fi
AC_SUBST([libcurl])

gcs=disabled
if test "$enable_gcs" != no; then
  if test $libcurl = enabled; then
    AC_DEFINE([ENABLE_GCS], 1, [Define if HTSlib should enable GCS support.])
    gcs=enabled
  else
    case "$enable_gcs" in
      check) AC_MSG_WARN([GCS support not enabled: requires libcurl support]) ;;
      *) AC_MSG_ERROR([GCS support not enabled

Support for Google Cloud Storage URLs requires libcurl support to be enabled
in HTSlib.  Configure with --enable-libcurl in order to use GCS URLs.])
      ;;
    esac
  fi
fi
AC_SUBST([gcs])

s3=disabled
if test "$enable_s3" != no; then
  if test $libcurl = enabled; then
    s3=enabled
    need_crypto="$enable_s3"
  else
    case "$enable_s3" in
      check) AC_MSG_WARN([S3 support not enabled: requires libcurl support]) ;;
      *) AC_MSG_ERROR([S3 support not enabled

Support for Amazon AWS S3 URLs requires libcurl support to be enabled
in HTSlib.  Configure with --enable-libcurl in order to use S3 URLs.])
      ;;
    esac
  fi
fi

CRYPTO_LIBS=
if test $need_crypto != no; then
  AC_CHECK_FUNC([CCHmac],
    [AC_DEFINE([HAVE_COMMONCRYPTO], 1,
               [Define if you have the Common Crypto library.])],
    [save_LIBS=$LIBS
     AC_SEARCH_LIBS([HMAC], [crypto],
       [AC_DEFINE([HAVE_HMAC], 1, [Define if you have libcrypto-style HMAC().])
        case "$ac_cv_search_HMAC" in
          -l*) CRYPTO_LIBS=$ac_cv_search_HMAC ;;
        esac],
     [case "$need_crypto" in
     check) AC_MSG_WARN([S3 support not enabled: requires SSL development files])
         s3=disabled ;;
     *) AC_MSG_ERROR([SSL development files not found

Support for AWS S3 URLs requires routines from an SSL library.  Building
HTSlib with libcurl enabled requires SSL development files to be installed
on the build machine; you may need to ensure a package such as libgnutls-dev,
libnss3-dev, or libssl-dev (on Debian or Ubuntu Linux, corresponding to the
libcurl4-*-dev package installed), or openssl-devel (on RPM-based Linux
distributions or Cygwin) is installed.

Either configure with --disable-s3 or resolve this error to build HTSlib.]) ;;
       esac])
     LIBS=$save_LIBS])
dnl Only need to add to static_LIBS if not building as a plugin
  if test "$enable_plugins" != yes ; then
     static_LIBS="$static_LIBS $CRYPTO_LIBS"
  fi
fi

dnl Look for PTHREAD_MUTEX_RECURSIVE.
dnl This is normally in pthread.h except on some broken glibc implementations.
AC_CHECK_DECL(PTHREAD_MUTEX_RECURSIVE, [], [AC_DEFINE([_XOPEN_SOURCE],[600], [Needed for PTHREAD_MUTEX_RECURSIVE])], [[#include <pthread.h>]])

if test "$s3" = enabled ; then
   AC_DEFINE([ENABLE_S3], 1, [Define if HTSlib should enable S3 support.])
fi
AC_SUBST([s3])
AC_SUBST([CRYPTO_LIBS])

AC_SUBST([pc_requires])
AC_SUBST([private_LIBS])
AC_SUBST([static_LDFLAGS])
AC_SUBST([static_LIBS])

AC_CONFIG_FILES([config.mk htslib.pc.tmp:htslib.pc.in])
AC_OUTPUT
