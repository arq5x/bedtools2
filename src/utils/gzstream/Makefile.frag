# ============================================================================
# gzstream, C++ iostream classes wrapping the zlib compression library.
# Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# ============================================================================
#
# File          : Makefile
# Revision      : $Revision: 1.3 $
# Revision_date : $Date: 2001/10/04 15:09:28 $
# Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
#
# ============================================================================

# ----------------------------------------------------------------------------
# adapt these settings to your need:
# add '-DGZSTREAM_NAMESPACE=name' to CPPFLAGS to place the classes
# in its own namespace. Note, this macro needs to be set while creating
# the library as well while compiling applications based on it.
# As an alternative, gzstream.C and gzstream.h can be edited.
# ----------------------------------------------------------------------------

BUILT_OBJECTS += obj/gzstream.o

# Uses $(CPPFLAGS) rather than $(ALL_CPPFLAGS) to avoid bedtools includes etc.

obj/gzstream.o: src/utils/gzstream/gzstream.C src/utils/gzstream/gzstream.h
	$(CXX) $(ALL_CXXFLAGS) -I$(<D) $(CPPFLAGS) -c -o $@ $<
