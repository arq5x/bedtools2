dnl     HOW TO CUT A NEW RELEASE
dnl
dnl  1) You will need (at least):
dnl     * Autoconf 2.69
dnl     * Automake 1.15
dnl     * Autoconf-archive 2016.09.16
dnl
dnl  2) Make your final commit
dnl  3) Tag the commit
dnl  4) Bootstrap the autotool toolchain, by running
dnl
dnl       autoreconf -vif
dnl
dnl  5) Run the ./configure script
dnl  6) Run 'make distcheck' in order to build a tarball and test it
dnl  7) Upload the resulting tarball to GitHub

AC_PREREQ([2.69])
AC_INIT(
	[bedtools],
	m4_esyscmd([./git-version-gen .tarball-version]),
	[aaronquinlan@gmail.com],
	[bedtools],
	[https://bedtools.readthedocs.io/])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src/bedtools.cpp])
AC_LANG([C++])



dnl ========================================================================
dnl Check whether we want to set defaults for CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ========================================================================
AC_MSG_CHECKING([whether configure should try to set CXXFLAGS/CPPFLAGS/LDFLAGS])
AS_IF([test "x${CXXFLAGS+set}" = "xset" || test "x${CPPFLAGS+set}" = "xset" || test "x${LDFLAGS+set}" = "xset"],
	[enable_flags_setting=no],
	[enable_flags_setting=yes]
)
AC_MSG_RESULT([${enable_flags_setting}])

AX_CHECK_ENABLE_DEBUG
AS_IF([test "x${enable_debug}" = "xno"],
	[
		cxxflags_test="-O2"
		cppflags_test=""

		AC_CANONICAL_HOST
		AS_CASE([${host_os}],
			[darwin*],
			[ldflags_test="-Wl,-dead_strip_dylibs"],
			[linux*],
			[ldflags_test="-Wl,-O1 -Wl,--as-needed -Wl,--sort-common -Wl,--gc-sections"]
		)
	]
)



dnl =======================================
dnl Check for standard headers and programs
dnl =======================================
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

AC_PROG_CXX
AC_PROG_SED
AC_PROG_RANLIB
AM_PROG_AR

AC_SYS_LARGEFILE

dnl Enable/disable bundled bamtools
AC_ARG_ENABLE([bundled-bamtools],
	AS_HELP_STRING([--disable-bundled-bamtools], [Disable bundled bamtools building; use system bamtools instead]))
AC_ARG_VAR([BAMTOOLS_CFLAGS], [Bamtools CFLAGS, used for passing in -Iincludedir when using --disable-bundled-bamtools])
AC_ARG_VAR([BAMTOOLS_LIBS], [Bamtools LIBS, used for passing in -Llibdir when using --disable-bundled-bamtools])

AM_CONDITIONAL([BUNDLED_BAMTOOLS], [test "x$enable_bundled_bamtools" != "xno"])

dnl Find zlib
PKG_CHECK_MODULES([ZLIB], [zlib])

dnl Documentation
AC_ARG_ENABLE([doc],
	AS_HELP_STRING([--enable-doc], [Enable building documentation (requires Sphinx)]))

AM_CONDITIONAL([ENABLE_DOCUMENTATION], [test "x$enable_doc" = "xyes"])

AS_IF([test "x$enable_doc" = "xyes"], [
	dnl Try to find sphinx-build
	AC_CHECK_PROGS([SPHINXBUILD], [sphinx-build sphinx-build3 sphinx-build2], [no])
	AS_IF([test "x$SPHINXBUILD" = "xno"], [
		dnl Die if sphinx-build cannot be found but documentation building was requested
		AC_MSG_ERROR([sphinx-build is required to build documentation])
	])
])



dnl ==================================
dnl Set CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ==================================
AS_IF([test "x${enable_flags_setting}" = "xyes" && test "x${enable_debug}" = "xno"],
	[
	dnl Optional flags that can be added: -Wall -Wextra
	AX_APPEND_COMPILE_FLAGS([${cxxflags_test}], [CXXFLAGS])
	CXXFLAGS="$( echo ${CXXFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g' )"

	AX_APPEND_COMPILE_FLAGS([${cppflags_test}], [CPPFLAGS])
	CPPFLAGS="$( echo ${CPPFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g' )"

	AX_APPEND_LINK_FLAGS([${ldflags_test}], [LDFLAGS])
	LDFLAGS="$( echo ${LDFLAGS}  | $SED -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g' )"
	]
)



dnl ===================
dnl Initialise Automake
dnl ===================
AM_INIT_AUTOMAKE([1.15 foreign dist-bzip2 no-dist-gzip subdir-objects silent-rules])
AM_SILENT_RULES([yes])



dnl ========
dnl Finalise
dnl ========
AC_CONFIG_FILES([
	Makefile
	docs/Makefile
])
AC_OUTPUT



dnl ============================
dnl Report configuration to user
dnl ============================
AC_MSG_RESULT([
	${PACKAGE_NAME} ${VERSION}
	Using ${ax_cv_cxx_compiler_vendor}, ${ax_cv_cxx_compiler_version}

	CXX:        ${CXX}
	CXXFLAGS:   ${CXXFLAGS}
	CPPFLAGS:   ${CPPFLAGS}
	LDFLAGS:    ${LDFLAGS}
])
